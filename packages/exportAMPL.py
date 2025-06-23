from amplpy import AMPL, Environment
from amplpy import Set
import pandas as pd

# Standard package imports
import pickle
import sys
import os
import re
import csv
from datetime import datetime
import gzip

# Interactive widget imports
import ipywidgets as widgets
from ipywidgets import interact, Layout
import pandas as pd
import numpy as np

# Plotting imports
from functools import partial
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Plotly theming
import plotly.io as pio

import sys
import heuristics as new_module



PAPER_LAYOUT=dict(font=dict(family="Computer Modern", size=18),
                  margin=dict(t=40))
pio.templates["paper"] = go.layout.Template(layout=PAPER_LAYOUT)
pio.templates.default = "paper"

if os.name == 'nt':
    # Windows, update with directory containing ampl.exe
    amplPath = Environment(os.path.normpath('C:/ampl'))

else:
    # MacOs, Linux typically just works
    amplPath = None


def padTruncate(input, target_len):
    """
    Extend or shorten the input list to the target length

    :param input: A list to be padded or truncated
    :param target_len: The ending length of the input list
    :return: List of target length
    """
    # Value to insert if padding needed
    pad_value = '_'
    
    return input[:target_len] + [pad_value]*(target_len - len(input))

# def getUniqueColumns(columns):
    # """Helper function to force unique dataframe column names"""
    # if not len(columns) == len(set(columns)):
        # columns = [col+f"({i})" for i,col in enumerate(columns)]
        
    # return columns

def transposeTupleList(tupList):
    """
    Transpose a row-wise list of tuples into a column-wise list of tuples.
    Example:
        transposeTupleList([(1, 3), (2, 4)])
        returns [(1, 2), (3, 4)]

    :param tupList: A list of tuples
    :return: A transposed list of tuples
    """
    return list(zip(*tupList))
    
def interact_row_wrap(_InteractFactory__interact_f=None, **kwargs):
    def patch(obj):
        if hasattr(obj.widget, 'layout'):
            obj.widget.layout = Layout(flex_flow='row wrap')
        return obj
    if _InteractFactory__interact_f is None:
        def decorator(f):
            obj = interact(f, **kwargs)
            return patch(obj)
        return decorator
    else:
        obj = interact(_InteractFactory__interact_f, **kwargs)
        return patch(obj)
        

def scaler(s):
    if s.max() > s.min():
        return (s - s.min()) / (s.max() - s.min())
    else:
        a = s.copy()
        a[:] = 1
        return a
    
def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)
    
def get_options(s):
    # create {label: (label, name, slice, suffix)} entries for all_names
    options = dict()
    
    for name in s.all_names:
        if s.data[name] is None:
            continue
        count = s.info[name]['indexarity']
        offset = -1 if (name not in s.sets) and (s.data[name].index.nlevels < 2) else 0
        
        if count == 0:
            if name in s.sets:
                options[name] = (name, name, -1)
            else:
                if name in s.cons:
                    options[name+'.lhs'] = (name+'.lhs', name, -2, 'lhs')
                    options[name+'.rhs'] = (name+'.rhs', name, -2, 'rhs')
                else:
                    options[name] = (name, name, -2)
        else:
            for i in range(count):
                temp = ['*'] * count
                temp[i] = '_'
                outStr = ','.join(temp)
                if name in s.cons:
                    options[f"{name}[{outStr}].lhs"] = (f"{name}[{outStr}].lhs", name, i+offset, 'lhs')
                    options[f"{name}[{outStr}].rhs"] = (f"{name}[{outStr}].rhs", name, i+offset, 'rhs')
                else:
                    options[f"{name}[{outStr}]"] = (f"{name}[{outStr}]", name, i+offset)
    return options

def sliceSeries(s, name, level, suffix=None):
    if level < 0:
        if suffix is None:
            return [(0, s.data[name]),]
        return [(0, s.data[name][suffix]),]
    if suffix is None:
        return [(i, s.data[name].xs(i, level=level)) for i in s.data[name].index.levels[level].to_list()]
    else:
        return [(i, s.data[name][suffix].xs(i, level=level)) for i in s.data[name].index.levels[level].to_list()]

def getAMPLEntity(entityList, extraMethods=None, suffixes=None, names=None):
    """
    Collect all AMPL entities in entityList to a dictionary.

    :param entityList: A list of name, value pairs where name is the name of an AMPL
        entity, and value is a handle to the AMPL object
    :param extraMethods: A list of methods to call on the entity (in addition to the defaults)
    :param suffixes: A list of .suffixes to add to the getValues() call on the entity
    :return: A dictionary whose keys are the names of the entities, and values are a dictionary 
        of the outputs of the methods specified plus data
    """
    # create an empty dict
    entities = {}
    
    # default methods
    methods = ['getIndexingSets', 'indexarity', 'isScalar', 'name', 'numInstances', 'toString'] 
    
    # add extra methods (repeats are included once)
    if extraMethods is not None:
        methods = list(set(methods + extraMethods))
    
    # for each of the entities (as in for each parameter, variable, etc.)
    for entity in entityList:
        try:
            # split the tuple returned by ampl
            name, val = entity
            
            if names is not None and name not in names:
                continue
                
            entities[name] = {}
            
            # call each of the methods provided as input (gets various attributes of the entity)
            for method in methods:
                entities[name][method] = getattr(val, method)()
                
            try:
                # Scalar parameters, variables, sets, and objectives
                if entities[name]['isScalar']:
                    if suffixes is None:
                        data = getattr(val, 'getValues')().toList()
                        temp = pd.Series(data)
                        
                    else:
                        # likely only used for constraints
                        data = {key:val for key,val in zip(suffixes,
                                                           transposeTupleList(getattr(val, 'getValues')(suffixes).toList()))}
                        temp = pd.DataFrame(data)
                        
                    temp.index += 1
                    entities[name]['data'] = temp
                
                # non-scalar parameters, variables, and sets
                else:
                    # for indexed sets, which are different than everything else
                    if isinstance(val, Set):
                        idx_list = []
                        val_list = []

                        # for each unique index of the set, get the index and the set values
                        for subSet in list(val):
                            if isinstance(subSet[0], tuple):
                                idx_list.append(subSet[0])
                            else:
                                idx_list.append(subSet[:-1])
                            temp = pd.Series(subSet[-1].getValues().toList())
                            temp.index += 1
                            val_list.append(temp)
                        
                        # assemble into a dataframe with columns for each index and one of the value
                        idx = pd.MultiIndex.from_tuples(idx_list)
                        idx.names = padTruncate(entities[name]['getIndexingSets'], idx.nlevels)
                        entities[name]['data'] = pd.Series(val_list, index=idx, name=name)

                    else:
                        if suffixes is None:
                            # for non-scalar parameters and variables, assemble the dataframe
                            data_list = transposeTupleList(val.getValues().toList())
                            idx = pd.MultiIndex.from_tuples(transposeTupleList(data_list[:-1]))
                            idx.names = padTruncate(entities[name]['getIndexingSets'], idx.nlevels)
                            entities[name]['data'] = pd.Series(data_list[-1], index=idx, name=name)
                            
                        else:
                            # for non-scalar constraints 
                            data_list = transposeTupleList(val.getValues(suffixes).toList())
                            idx = pd.MultiIndex.from_tuples(transposeTupleList(data_list[:-len(suffixes)]))
                            idx.names = padTruncate(entities[name]['getIndexingSets'], idx.nlevels)
                            data = {key:val for key,val in zip(suffixes, data_list[-len(suffixes):])}
                            entities[name]['data'] = pd.DataFrame(data, index=idx)

            # throws a TypeError if no value is found (e.g., variables may not have values if the model is not solved)
            except TypeError as e:
                print('## Get entity error ##')
                print(name, e)
                # continue
                # raise
        except RuntimeError as e:
            # just go to the next one if an error occurs
            print('## Get entity error ##')
            print(name, e)
            # raise
        
        # many of the amplpy methods do not reasonably handle even trivial errors
        #   this was added in testing, it is unclear if it is still needed
        except Exception as e:
            # just go to the next one if an error occurs
            print('## Get entity error ##')
            print(name, e)
            # raise
    
    # return the collected dict
    return entities
    

def parseConstraints(conDict, conList, ampl, names=None):
    """
    Calculate the left and right hand side constraint values (useful for plotting and debugging)
    
    :param conDict: A dictionary of dictionaries describing a model's constraints (output of getAMPLEntity)
    :param conList: A list of name, value pairs where name is the name of an AMPL
        entity, and value is a handle to the AMPL object
    :param ampl: An instance of amplpy.AMPL, which is an active connection to ampl
    :return conDict: A updated nested dictionary with additional information and columns in 'data'
    """
    # walk through the constraints and parse the constraint equations
    for name, constraint in conList:
        
        if names is not None and name not in names:
            continue
    
        # get the constraint as a string
        conStr = conDict[name]['toString']
        
        # get the indexing set as a string
        temp = re.search(r'{(.|\s)*?}', conStr)
    
        if temp is not None and len(conDict[name]['getIndexingSets']):
            idxSet = temp.group(0)
        else:
            idxSet = ''
        
        # find the colon at the beginning of the constraint equation
        endDefStr = re.search(r'(:)(?![^{]*\})', conStr)
        b, e = endDefStr.span() # beginning and end of everything before the equation
        
        # search for the inequality operator in the constraint equation
        operStr = re.search(r'(<=|==|>=|=)(?!([^{(i])*([\}\)]|(\sthen\s)))', conStr[e:])
        ob, oe = operStr.span()
        ob, oe = ob+e, oe+e # Beginning and end of the operator
        
        # left and right hand side (lhs and rhs) of the constraint equation
        lhs = conStr[e:ob].strip().replace('\n','')
        rhs = conStr[oe:-1].strip().replace('\n','')
        
        # strings showing how the constraint was parsed
        conDict[name]['lhs'] = lhs + conStr[ob:oe] # Add the operator to the end of the LHS
        conDict[name]['rhs'] = rhs
        conDict[name]['idx_set'] = idxSet
        
        if 'data' not in conDict[name].keys():
            continue
        
        # ask ampl to evaluate the equations to values, return a panda dataframe
        try:
            data_list = ampl.getData(idxSet+lhs, idxSet+rhs).toList()
            data = transposeTupleList(data_list)
            
            if conDict[name]['isScalar']:
                offset = 0
            else:
                offset = 1

            conDict[name]['data']['lhs'] = data[offset]
            conDict[name]['data']['rhs'] = data[offset+1]
        
        # throws a KeyError if no value found
        except KeyError as e:
            print('## Constraint parsing error ##')
            print(name, e)
            
        except TypeError as e:
            print('## Constraint parsing error ##')
            print(name, e)
            
        except Exception as e:
            print('## Constraint parsing error ##')
            print(name, e)
        
    return conDict
    
    
def writeConsReport(Files, Cons):
    """
    Write out a csv file showing how model constraints were parsed (the regular expression used is not extensively tested)
    
    :param Files: A dictionary of solve meta data (files used, solver log, time, etc.)
    :param Cons: A dictionary of dictionaries describing a model's constraints (output of parseConstraints)
    :return: None
    """
    # report file name (for checking constraint parsing)
    reportFile = 'constraintReport.csv'
    
    if (sys.version_info > (3, 0)):
        # python 3 convention
        kwargs = {'mode': 'w', 'newline': ''}
        
    else:
        # python 2 convention
        kwargs = {'mode': 'w'}

    # open the file for writing
    try:
        with open(reportFile, **kwargs) as f:
            # csv writer object
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            
            if Files is not None:
                # write the system time, and ampl file names
                writer.writerow(["'" + Files['time']])
                writer.writerow([Files[key] for key in Files.keys() if key != 'time'])
            writer.writerow([])
            
            # write header line
            writer.writerow(['Name', 'Full Constraint', 'Indexing Set', 'LHS', 'RHS'])
            
            # for each named constraint
            for constraint in sorted(list(Cons.keys())):
                # index the items we parsed earlier
                conStr = Cons[constraint]['toString']
                idxSet = Cons[constraint]['idx_set']
                lhs = Cons[constraint]['lhs']
                rhs = Cons[constraint]['rhs']
                
                # write them out
                writer.writerow([Cons[constraint]['name'], conStr, idxSet, lhs, rhs])
                
        print('Constraint report written to ' + reportFile)
                
    except Exception as e:
        print(e)
        print('WARNING cannot access ' + reportFile)


class amplSolution():
    """
    An object to hold all information from an AMPL solution.
    """
    
    def __init__(self, Sets, Params, Vars, Objs, Cons, Files=None):
        """
        :param Sets: A dictionary of dictionaries describing AMPL sets (output of getAMPLEntity)
        :param Params: A dictionary of dictionaries describing AMPL parameters (output of getAMPLEntity)
        :param Vars: A dictionary of dictionaries describing AMPL variables (output of getAMPLEntity)
        :param Objs: A dictionary of dictionaries describing AMPL objective(s) (output of getAMPLEntity)
        :param Cons: A dictionary of dictionaries describing AMPL constraints (output of parseConstraints)
        :param Files: A dictionary of solve meta data (files used, solver log, time, etc.)
        """
        # add current time if no meta-data is given
        if Files is None:
            Files = {'time': str(datetime.now())}
        
        # store meta_daa
        self.meta_data = Files
        
        # store names of AMPL entities
        self.sets = list(Sets.keys())
        self.params = list(Params.keys())
        self.vars = list(Vars.keys())
        self.objs = list(Objs.keys())
        self.cons = list(Cons.keys())
        
        # all entity names
        self.all_names = [item for sublist in [Sets, Params, Vars, Objs, Cons]
                               for item in sublist]
        
        # all info dictionaries
        self.info = dict()
        [self.info.update(d) for d in [Sets, Params, Vars, Objs, Cons]]
        
        # move the 'data' objects to a different dictionary
        self.data = {key:self.info[key].pop('data', None) for key in self.all_names}

        # add convenience attributes for all the named AMPL entities
        self.update_attrs()
        
    def update_attrs(self):
        """
        Add convenience attributes for all the named AMPL entities, unless they already exist
            
        """
        for name in self.all_names:
            i = 1
            temp = name
            while temp in vars(self):
                temp = name + str(i)
                i += 1
                
            setattr(self, temp, self.data[name])
        
    def display(self, nameList):
        if not isinstance(nameList, list):
            nameList = [nameList]
            
        for name in nameList:
            print(f"{name}{self.info[name]['getIndexingSets']}")
            if self.data[name].index.nlevels > 1:
                print(self.data[name].unstack())
            else:
                print(self.data[name])
            print()
            
    def displayInfo(self, nameList, keys=None):
        if not isinstance(nameList, list):
            nameList = [nameList]
            
        for name in nameList:
            print(f"{name}")
            for key in self.info[name].keys():
                if keys is not None and key not in keys:
                    continue
                print(f"{key}: {self.info[name][key]}")
            print()
            
    
    def multiPlot(self, x=None, y=None, plot_type=None):
        options = {'None': None}
        options.update(get_options(self))

        plot_funcs = dict(scatter = partial(go.Scatter, mode='markers'),
                          line = partial(go.Scatter, mode='markers+lines'),
                          bar = go.Bar,
        #                   box = go.Box,
                         )
        plot_types = list(plot_funcs.keys())

        x_default = options[x] if x is not None else list(options.values())[0]
        y_default = tuple([options[yi] for yi in y]) if y is not None else (list(options.values())[1],)
        plot_default = plot_types[0]
        index_options = {'X and Y share an index': 0,
                         'Use X values to index Y[*]': 1,
                         'Ignore indexes, just plot': 2}

        switches = dict(
            x_og = widgets.Dropdown(options=options, value=x_default, description='X'),
            y_og = widgets.SelectMultiple(options=options, value=y_default, description='Y', rows=8),
            plot_type = widgets.Dropdown(options=plot_types, value=plot_default, description='Type'),
            index = widgets.Dropdown(options=index_options, value=0, description='Index'),
            scale = widgets.Checkbox(value=False, description='Scale Y'),
        )

        def on_change(change):
            if change['type'] == 'change' and change['name'] == 'value':
                if change['new'] is not None:
                    x_name = change['new'][1]
                    if x_name in self.sets:
                        switches['index'].value = 1
                    else:
                        switches['index'].value = 0

        switches['x_og'].observe(on_change)

        def amplPlot(x_og, index, plot_type, y_og, scale):    
            try:
                plot_func = plot_funcs[plot_type]
                fig = go.Figure()
                
                for y in y_og:
                    if x_og is None and y is None:
                        return fig
                    if x_og is None:
                        y_name = y[1]
                        x_label, x_name, x_count, x_data = ','.join(self.info[y_name]['getIndexingSets']), None, -2, [None]
                    else:
                        x_label, x_name, x_count = x_og[0], x_og[1], x_og[2]
                        x_data = sliceSeries(self, *x_og[1:])
               
                    if y is None:
                        y_label, y_name, y_count, y_data = ','.join(self.info[x_name]['getIndexingSets']), None, -2, [None]
                    else:
                        y_label, y_name, y_count = y[0], y[1], y[2]
                        y_data = sliceSeries(self, *y[1:])


                    if x_count < 0: x_data *= len(y_data)
                    if y_count < 0: y_data *= len(x_data)

                    for (i, xdat), (j, ydat) in zip(enumerate(x_data), enumerate(y_data)):
                        if xdat is None: 
                            index = 2
                            if isinstance(ydat[1].index, pd.MultiIndex):
                                idx = ydat[1].index.levels[0]
                            else:
                                idx = ydat[1].index
                            xdat = (0, idx)
                        if ydat is None: 
                            index = 2
                            if isinstance(xdat[1].index, pd.MultiIndex):
                                idx = xdat[1].index.levels[0]
                            else:
                                idx = xdat[1].index
                            ydat = (0, idx)

                        xidx, x = xdat
                        yidx, y = ydat

                        if x_name in self.sets and x_count >= 0: x = x.iloc[0] # hack for indexed sets
                        if y_name in self.sets and y_count >= 0: y = y.iloc[0] # hack for indexed sets

                        if y_count >= 0:
                            lbl = rreplace(y_label, '_', str(yidx), 1)
                        else:
                            lbl = y_label

                        if index == 2:
                            # ignore indexes and just plot
                            if scale: y = scaler(y)
                            fig.add_trace(plot_func(x=x, y=y, name=lbl))

                        elif index == 1:
                            # use values of X to index Y
                            if isinstance(y.index, pd.MultiIndex):
                                idx = pd.MultiIndex.from_arrays([x])
                            else: idx = x
                            y = y.reindex(idx)
                            if y_count == -1 and y_name not in self.sets:
                                lbl = rreplace(lbl, '_', rreplace(x_label, '_', str(xidx), 1), 1)
                            else: lbl = rreplace(lbl, '*', x_name, 1)
                            if scale: y = scaler(y)
                            fig.add_trace(plot_func(x=x, y=y, name=lbl))

                        else:
                            # X and Y are indexed by the same set
                            if isinstance(y.index, pd.MultiIndex):
                                if isinstance(x.index, pd.MultiIndex):
                                    idx = pd.MultiIndex.from_arrays([*x.index.levels])
                                else: idx = pd.MultiIndex.from_arrays([x.index.values])
                            else: 
                                if isinstance(x.index, pd.MultiIndex):
                                    idx = x.index.levels[0]
                                else: idx = x.index
                            y = y.reindex(idx)
                            if y_count == -1 and y_name not in self.sets and y_name not in self.cons:
                                lbl = rreplace(lbl, '_', str(xidx), 1)
                            if scale: y = scaler(y)
                            fig.add_trace(plot_func(x=x, y=y, name=lbl))
                            
                fig.update_xaxes(title=x_label)
                fig.update_yaxes(title=y_label)
                fig.update_layout(width=900, hovermode='x', barmode='group')
                return fig

            except:
                print('Something broke...')
                raise
                
        return interact_row_wrap(amplPlot, **switches)
    
            
    def plot(self, x=None, y=None, plot_type=None):
        options = {'None': None}
        options.update(get_options(self))

        plot_funcs = dict(scatter = partial(go.Scatter, mode='markers'),
                          line = partial(go.Scatter, mode='markers+lines'),
                          bar = go.Bar,
        #                   box = go.Box,
                         )
        plot_types = list(plot_funcs.keys())

        x_default = options[x] if x is not None else list(options.values())[0]
        y_default = options[y] if y is not None else list(options.values())[1]
        plot_default = plot_type if plot_type is not None else plot_types[0]
        index_options = {'X and Y share an index': 0,
                         'Use X values to index Y[*]': 1,
                         'Ignore indexes, just plot': 2}

        switches = dict(
            x = widgets.Dropdown(options=options, value=x_default, description='X'),
            y = widgets.Dropdown(options=options, value=y_default, description='Y'),
            plot_type = widgets.Dropdown(options=plot_types, value=plot_default, description='Type'),
            index = widgets.Dropdown(options=index_options, value=0, description='Index')
        )

        def on_change(change):
            if change['type'] == 'change' and change['name'] == 'value':
                if change['new'] is not None:
                    x_name = change['new'][1]
                    if x_name in self.sets:
                        switches['index'].value = 1
                    else:
                        switches['index'].value = 0

        switches['x'].observe(on_change)

        def amplPlot(x, y, plot_type, index):               
            try:
                plot_func = plot_funcs[plot_type]
                fig = go.Figure()
                
                if x is None and y is None:
                    return fig
                if x is None:
                    y_name = y[1]
                    x_label, x_name, x_count, x_data = ','.join(self.info[y_name]['getIndexingSets']), None, -2, [None]
                else:
                    x_label, x_name, x_count = x[0], x[1], x[2]
                    x_data = sliceSeries(self, *x[1:])
                   
                if y is None:
                    y_label, y_name, y_count, y_data = ','.join(self.info[x_name]['getIndexingSets']), None, -2, [None]
                else:
                    y_label, y_name, y_count = y[0], y[1], y[2]
                    y_data = sliceSeries(self, *y[1:])
                
                if x_count < 0: x_data *= len(y_data)
                if y_count < 0: y_data *= len(x_data)
                    
                for (i, xdat), (j, ydat) in zip(enumerate(x_data), enumerate(y_data)):
                    if xdat is None: 
                        index = 2
                        if isinstance(ydat[1].index, pd.MultiIndex):
                            idx = ydat[1].index.levels[0]
                        else:
                            idx = ydat[1].index
                        xdat = (0, idx)
                    if ydat is None: 
                        index = 2
                        if isinstance(xdat[1].index, pd.MultiIndex):
                            idx = xdat[1].index.levels[0]
                        else:
                            idx = xdat[1].index
                        ydat = (0, idx)
                        
                    xidx, x = xdat
                    yidx, y = ydat
                        
                    if x_name in self.sets and x_count >= 0: x = x.iloc[0] # hack for indexed sets
                    if y_name in self.sets and y_count >= 0: y = y.iloc[0] # hack for indexed sets
                    
                    if y_count >= 0:
                        lbl = rreplace(y_label, '_', str(yidx), 1)
                    else:
                        lbl = y_label
                        
                    if index == 2:
                        # ignore indexes and just plot
                        fig.add_trace(plot_func(x=x, y=y, name=lbl))
                        
                    elif index == 1:
                        # use values of X to index Y
                        if isinstance(y.index, pd.MultiIndex):
                            idx = pd.MultiIndex.from_arrays([x])
                        else: idx = x
                        y = y.reindex(idx)
                        if y_count == -1 and y_name not in self.sets:
                            lbl = rreplace(lbl, '_', rreplace(x_label, '_', str(xidx), 1), 1)
                        else: lbl = rreplace(lbl, '*', x_name, 1)
                        fig.add_trace(plot_func(x=x, y=y, name=lbl))

                    else:
                        # X and Y are indexed by the same set
                        if isinstance(y.index, pd.MultiIndex):
                            if isinstance(x.index, pd.MultiIndex):
                                idx = pd.MultiIndex.from_arrays([*x.index.levels])
                            else: idx = pd.MultiIndex.from_arrays([x.index.values])
                        else: 
                            if isinstance(x.index, pd.MultiIndex):
                                idx = x.index.levels[0]
                            else: idx = x.index
                        y = y.reindex(idx)
                        if y_count == -1 and y_name not in self.sets and y_name not in self.cons:
                            lbl = rreplace(lbl, '_', str(xidx), 1)
                        fig.add_trace(plot_func(x=x, y=y, name=lbl))
                            
                fig.update_xaxes(title=x_label)
                fig.update_yaxes(title=y_label)
                fig.update_layout(width=900, hovermode='x', barmode='group')
                return fig
                
            except:
                print('Something broke...')
                # raise
                    
        return interact_row_wrap(amplPlot, **switches).widget   
            

def exportAll(ampl, Files=None, writeReport=False, names=None):
    """
    Export all model components to nested dictionary variables
    
    :param ampl: An instance of amplpy.AMPL, which is an active connection to ampl
    :param Files: An optional dictionary containing meta data
    :param writeReport: Boolean if the constraint report should be written
    :return amplSolution: An object containing all parts of the solution to an AMPL model
    """
    try:
        # get sets
        setMethods = ['arity']
        setList = list(ampl.getSets())
        Sets = getAMPLEntity(setList, extraMethods=setMethods, names=names)
        
        # get parameters
        paramMethods = ['hasDefault', 'isSymbolic']
        paramList = list(ampl.getParameters())
        Params = getAMPLEntity(paramList, extraMethods=paramMethods, names=names)
        
        # get variables
        varList = list(ampl.getVariables())
        Vars = getAMPLEntity(varList)
        
        # get objective(s)
        objMethods = ['astatus', 'exitcode', 'message', 'minimization', 'result', 'sstatus']
        objList = list(ampl.getObjectives())
        Objs = getAMPLEntity(objList, extraMethods=objMethods, names=names)
        
        # get constraints
        conList = list(ampl.getConstraints())
        suffixes = ['lb', 'body', 'ub', 'status', 'dual']
        Cons = getAMPLEntity(conList, suffixes=suffixes, names=names)
        Cons = parseConstraints(Cons, conList, ampl, names=names)
  
    # raise any exception (many things can go wrong)
    except Exception as e:
        # print(e)
        raise
        
    # if writeReport:
    #     # write out CSV to check constraint parsing
    #     writeConsReport(Files, Cons)
    
    
    return amplSolution(Sets, Params, Vars, Objs, Cons, Files)
    
def loadSolution1(filename):
    sys.modules['packages'] = new_module
    
    with gzip.open(filename, 'rb') as f:
        # read the solution variables
        solution = pickle.load(f)
        return solution
def loadSolution(filename):
    """
    Load a pickled solution from file
    
    :param filename: Name of a file containing a pickled solution 
    return 
    """
    
    # open the file for reading
    with gzip.open(filename, 'rb') as f:
        # read the solution variables
        solution = pickle.load(f)
    
    # add convenience attributes (if the loaded solution is an amplSolution object
    if isinstance(solution, amplSolution):
        solution.update_attrs()
    
    return solution

def saveSolution(ampl, fileName=None, meta_data=None, names=None):
    """
    Save an AMPL solution out to file from an AMPL connection
    
    :param ampl: An instance of amplpy.AMPL, which is an active connection to ampl
    :param Files: An optional dictionary containing meta data
    :return solution: an amplSolution object
    """
    # default values for the Files dict
    if meta_data is None:
        meta_data = {'pickleFile': 'solution.pkl',
                     'time': str(datetime.now())}
    
    # at least add the pickle file name
    elif isinstance(meta_data, dict) and 'pickleFile' not in meta_data.keys():
        meta_data['pickleFile'] = fileName
        
    solution = exportAll(ampl, Files=meta_data, writeReport=True, names=names)
    
    # open the file for writing and write out) protocol 2 used for compatibility from Linux to Windows
    #  and reasonable file size
    with gzip.open(fileName, 'wb') as f:
        pickle.dump(solution, f)
        
    return solution
    
    
def main(runFile, fileName):
    """
    Execute an AMPL run file 
    
    :param runFile: The name of an AMPL run file to be executed
    :param fileName: The name of a pickle file to write the solution out to 
    """
    # Open the connection to AMPL
    ampl = AMPL()
    
    meta_data = {'runFile': runFile, 
                 'pickleFile': fileName,
                 'time': str(datetime.now())}

    # Read the model file (and solve if solved in run file)
    ampl.read(runFile)
    
    saveSolution(ampl, fileName, meta_data)
    
    
"""If run at the command line, call main with command line arguments"""
if __name__ == '__main__':
    argc, argv = len(sys.argv), sys.argv
    
    # If not called with at least 3 input args, return and print the doc string above
    if argc < 2:
        print("Usage: python saveSolution.py <runFile> <opt pickleFile>")
        raise()


    # Read the run file name (assumes they are in the current directory)
    # Pickle file name (default is solution.pkl, or use name provided at command line)
    runFile, pickleFile = argv[1], argv[2] if argc == 3 else 'solution.pkl'
    
    main(runFile, pickleFile)
        
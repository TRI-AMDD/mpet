#!/usr/bin/env python
"""VKMLgui.py
""" 
 # ###################################################################
 #  FILE: "VKMLgui.py" 
 #
 #  Author: Alex Bartol <nanohub@alexbartol.com>
 #       Copyright 2009
 #
 # 
 # VKMLgui.py is a utility to allow for user input on the fly with a simple
 # understandable GUI or on the command line.
 #
 ## Variable types supported : integer, float, string, boolean, list, tuple, 
 ##                             file
 #
 # Module's graphical user mode can be used by using the --gui flag 
 # No flag is makes the program start run the callback with defaults and close
 #
 # when a parameter is created and given a default value, if the program was
 # run with that variable inputed on command line (--variable=12) then the 
 # default value given will be overriden by the commmand line argument
 # ###################################################################
 #


import inspect, time
from sys import argv as ARGS
from sys import stdout as stdout
from copy import deepcopy
from threading import Thread
import sys

from VKML.gui.tk.Parameter import Parameter
from VKML.gui.tk.Menu import Menu
from VKML.gui.tk.Gui import _GetParams



class Parser:
    """ Parser is designed to be an easy way to get input from the typical
    user. Instead of relying on the user being able to edit hard-coded
    variables or remembering verbose command line arguments, this class 
    generates a GUI on the fly for the user (or asks for variables
    on the command-line)
    """
    ON = 1
    OFF = -1
    AUTO = 0
    def __init__(self, title='',help_text='', interactive=ON, update_text='Update', **kwargs):
        """Starts a new Parser
        title is the desired title for the GUI window
        interactive modes are (ON, OFF, and AUTO)
        **kwargs - are potential callback functions
        """
        self.update_text = update_text
        self.params = []
        self.tabs = []
        self.gui_title = title
        self.interactive = interactive
        self.functions = []
        self.pre_functions = []
        self.post_functions = []
        self.ask_quit = True
        self.log = None
        self.log_on = False
        self.text = ''
        self.credits2 = ''
        self.help_text = help_text
        for arg in kwargs:
#            print arg
            if type(arg) == function:
                self.add_command(arg)

    def __call__(self):
        """More versitility
        """
        self.run()
        
    def set_credits(self, credits):
        """This is designed to allow the programmer to set the About menu item
         to include his credits in addition to mine
        """
        self.credits2 = credits
        
    def add(self, param):
        """ This function adds a variable to the parser.
        First argument must be a type of Parameter
        """
        self.params.append(param)

    def add_menu(self, menu):
        """This function adds an entire Menu into the parser
        First argument must be a type of Menu
        """
        self.tabs.append(menu)
        for param in menu.params:
            self.params.append(param)
        
    def add_pre_command(self, function):
        """Adds a command to be done FIRST and never to be repeated
        """
        kwargs = {}
        
        variables = inspect.getargspec(function)[0]
        
        for param in self.params:
            if param.name in variables:
                kwargs[param.name] = param
        self.pre_functions.append(self._command(function, kwargs, 
                                                self.interactive))

    def add_post_command(self, function):
        """Adds a command to be done LAST 
        (right before exit and not before then)
        """
        kwargs = {}
        
        variables = inspect.getargspec(function)[0]
        
        for param in self.params:
            if param.name in variables:
                kwargs[param.name] = param
        self.post_functions.append(self._command(function, kwargs,
                                 interactive=self.interactive))
        
    def run_post_commands(self):
        """runs all the commands designated to be done after the simulation
        takes place
        """
        for func in self.post_functions:
            func(forced=True)
    
    def run_pre_commands(self):
        """runs all of the commands designated to be done before the simulation
        takes place
        """
        for func in self.pre_functions:
            func(forced=True)
        
    def init_auto_mode(self):
        """checks for threashold timing (also runs the normal commands for the
        first time)
        """
        if self.interactive > self.OFF:
            for func in self.functions:
                func.init_call()
            
    def add_command(self, function):
        """function - The callback function
        
        This function determins which parameters are used and creates the 
        callback class accordingly
        """
        kwargs = {}
        
        #print str(len(self.params)) + " Number of args"
        variables = inspect.getargspec(function)[0]
        if inspect.getargspec(function)[2] == 'kwargs':
            for param in self.params:
                if param.type != 'label' and param.type != 'image':
                    kwargs[param.name] = param
                    #print param.name
        else:
            for param in self.params:
                if param.name in variables and param.type != 'label':
                    kwargs[param.name] = param
        
        self.functions.append(self._command(function, kwargs, self.interactive))
       

    def is_gui_mode(self):
        """Returns True/False if the mode is creating a gui or not
        """
        return '--gui' in sys.argv

    def is_text_mode(self):
        """Returns True/False if the mode is text based (no gui)
        """
        return not is_gui_mode(self)

    def get_mode(self):
        """Returns the string gui/text depending on mode
        """
        if self.is_gui_mode():
            return 'gui'
        return 'text'

    def set_mode(self, mode):
        """sets the mode if the mode is the string 'gui'/'text'
        """
        if mode.lower() == 'gui':
            sys.argv.append('--gui')
        if mode.lower() == 'text':
            while '--gui' in sys.argv:
                sys.argv.remove('--gui')

    def update(self, forced=False):
        """Calls the callback for all the functions available
        """
        for func in self.functions:
            func.update(forced=forced)
            
    def update_button(self):
        """This overrides the interactive mode because the button should force 
        an update
        """
        for func in self.functions:
            func()
            
    def __delete_duplicates__(self):
        """Steps through the parser and removes all duplicate entries
        """
        length = len(self.params)
        for i in range(length-1):
            for j in range(i+1, length):
                if i < length and j < length and self.params[i].name == self.params[j].name:
                    del self.params[j]
                    length-=1
    # def close_parser(self):
    #     (tk.Tk())

    
    def run(self):
        """ Starts the GUI, if gui is unavaliable or unwanted, it will start
        asking for commands on the command line instead
        """
        #determines parameters from command line arguments
        self.__delete_duplicates__()
        gui = False
        interactive = True
        command_line = False
        self.run_pre_commands()
        for arg in ARGS:
            if arg.lower() == '--gui':
                gui = True
            elif arg.lower() == '--quit':
                #self.ask_quit = False
                sys.exit()
            elif arg.lower() == '--help':
                print (self.help_text) 
                sys.exit()
#            elif arg.startswith('--help'):
#				if ARGS[0]=='dualfoil.py':
#					RunDualfoilHelp()
#					sys.exit()
#				elif ARGS[0]=='visualize.py':
#					RunVisualizeHelp()
#				sys.exit()
            elif arg.startswith('--') and '=' in arg:
                arg = arg[2:]
                name = arg.split('=')[0]
                value = arg.split('=')[-1]
                for param in self.params:
                    if name == param.name:
                        if param.type == float or param.type == int or param.type == str:
                            param.default = param.type(value)
                            param.out = param.default
                        if param.type == 'file':
                            param.filename = value
                            if param.filetype.startswith('o'):
                                param.out = open(param.filename, 'r')
                            elif param.filetype.startswith('s'):
                                param.out = open(param.filename, 'w')
                            else: #direcotry
                                param.out = param.filename
                        if type(param.type) == list:
                            def makeList(ls):
                                #print ls
                                if ls.count('[') == 1 and ls.count(']') == 1:
                                    ls = ls.replace('[', '').replace(']', '')
                                    rtn = ls.split(',')
                                    for i in range(len(rtn)):
                                        try:
                                            rtn[i] = int(rtn[i])
                                        except:
                                            try:
                                                rtn[i] = float(rtn[i])
                                            except:
                                                pass
                                    return rtn
                                bc = 0
                                rtn = []
                                for i in range(len(ls)):
                                    if ls[i] == '[':
                                        bc += 1
                                    if ls[i] == ']':
                                        bc -= 1
                                    if ls[i] == ',' and bc == 1:
                                       # print 'tick'
                                        ls = ls[:i] + '`' + ls[i+1:]
                                if bc != 0:
                                    print ('Incorrect format for type List')
                                    sys.exit()
                                ls = ls[1:]
                                ls = ls[:-1]
                                for i in ls.split('`'):
                                    rtn.append(makeList(i))
                                return rtn
                                
                            out = makeList(value)
                            param.out, param.default = out, out
        if gui:
            self.init_auto_mode()
            gui = _GetParams(self.tabs, parser=self, credits=self.credits2, title=self.gui_title, update_text=self.update_text)
            gui.get_result()
        else:
            #calls callback function before returning
            self.update(forced=True)
            return
            
    class _command:
        """Internal class designed to handle all callback functions
        """
        def __init__(self, func, kwargs, interactive):
            """func is the function that needs to be called
            kwargs is the dictionary that holds the the parameters
            interactive is the variable that states how often this should be updated 
            """
            self.f = func
            self.kwargs = kwargs
            self.last = {}
            self.interactive = interactive
            self.wid = -1
            self.time_threashold = 100
            self.last_time = 0
            self.asked = False
            self.out = None
            self.duration = 0
            for arg in self.kwargs.items():
                self.last[arg[0]] = arg[1]
           
        def update(self, forced=False):
            """Checks what mode it is currently in and runs functions 
            accordingly
            """
            if forced:
                self.force_run()
                return None
            if self.interactive == Parser.ON or \
                    (self.interactive == Parser.AUTO and \
                    self.duration < self.time_threashold):
                self()
		#self.start()
            elif self.duration >= self.time_threashold and not self.asked:
                self.asked = True
                self.interactive = False      
                    
        def set_last(self):
            """Records the previous values
            """
            for item in self.kwargs:
                if not self.kwargs[item].type == file:
                    #print self.last[item]()
                    self.last[item] = deepcopy(self.kwargs[item]())
                else: #FILE
                    self.last[item] = self.kwargs[item]()
                #print self.last[item]
                                
        def changed(self):
            """checks if there is a change, if change occurs returns True 
            """
            rtn = False
            for item in self.kwargs:
                if self.kwargs[item]() == self.last[item]:
                    continue
                else:
                    
                    self.set_last()
                    rtn = True
            self.set_last
            #self.update()
            return rtn
        
        def init_call(self): 
            """Calls the function (to start the gui in some cases), in order
            to establish a time that the simulation will take. 
            """
            self.changed() #needed to set self.last
            start_time = time.time()
            self.out = self.f(**self.last)
            self.duration = time.time() - start_time
            
        def run(self):
            self()
        
        def force_run(self):
            self(forced=True)
            
        def __call__(self, forced=False):
            """calls the function with arguments
            """
            if self.changed() or forced:
                self.out = self.f(**self.last)
                return self.out
                                                       
if __name__ == '__main__':
    """The following code is only a sample/test case for my software
    """
    P = Parser(title='Sample GUI', interactive=-1)
    
    t1 = Menu(title='tab1', parser=P)
    t2 = Menu(title='tab2', parser=P)
 
    def strfunc(**kwargs):
        for arg in kwargs:
            print (kwargs[arg])

    Parameter(name='string', menu=t1, default='text goes here', variable=str)
    Parameter(name='matrix1', menu=t2, default=[[1,2,3],[1,2,3],[1,2,3]], variable=list)
    Parameter(name='integer1', menu=t2, default=10, variable=int,interval=(-100, 100))
    
    Parameter(name='bool1', default=False, menu=t1, variable=bool)

    P.add_command(strfunc)
                                        
    P()
    
    

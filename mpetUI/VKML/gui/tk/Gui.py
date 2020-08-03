"""Gui.py
"""
 # ############################################################################
 # FILE: Gui.py
 #
 #
 #  Author: Alex Bartol <nanohub@alexbartol.com>
 #       Copyright 2009
 #
 # ############################################################################
try:
  #python 3
  import tkinter as tk
  from tkinter import messagebox as tkMessageBox
except ImportError:
  #python 2
  import Tkinter as tk
  import tkMessageBox
  
class _GetParams:
    """GetParams takes in a list of variables that it puts into GUI form
    The class only modifies the <Parameter>.out variable
    """
    def __init__(self, menus, parser, credits='', title='', update_text='Update'):
        """initializes variables and starts the GUI window
        
        Initialization Information:
        __init__(self, )
       
        menus - a list of type Menu that should be in the GUI
        parser - the parser that the program is using
        title - a string that will be the title bar for the GUI
        """
        self.results = []
        
        self.args = []
        self.update_text = update_text
        self.credits = credits
        for arg in self.args:
            self.results.append(arg.out)
        self.window = tk.Tk()
        self.window.withdraw() #this gets rid of the annoying TK window
        self.tabs = menus
        self.parser = parser
        self.entries = [] #the list of gui text boxes
        self.top = tk.Toplevel(self.window, padx=20, pady=20)
        self.top.columnconfigure(0, weight=1)
        self.top.rowconfigure(1, weight=1)
        #self.top.minsize(400,200)
        self.title = title
        if len(menus) == 1:
            self.top.title(self.title)
        else:
            try:
                self.top.title(self.title + ' - ' + menus[0].title) #sets title
       	    except IndexError:
       	        self.top.title(self.title)
        self.boolstatus = [] #keeps track of boolean values
        
        self._setup_menubar()
        self._setup_gui()
        
        
    def _setup_menubar(self):
        """Creates the Menubar and adds it to the GUI
        """
        self.menubar = tk.Menu(self.top)
        
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label='Quit', command=self.close)
        self.menubar.add_cascade(label='File', menu=filemenu)
        
        
        tabmenu = tk.Menu(self.menubar, tearoff=0)
        
        def set_value(tab):
            self.top.title(self.title + ' - ' + tab.title)
            for t in self.tabs:
                t.hide()
            i = tab.show(self.top)
            if i > 0:
                try:
                    self.u_button.destroy()
                except:
                    pass
                # self.u_button = tk.Button(self.top, text=self.update_text,
                #                           command=self._update_button)
                # self.u_button.grid(row=i+3, sticky=tk.S,
    	           #  column=0, columnspan=3)
                # self.u_button.bind('<Return>', self._update_button)
            else:
                try:
                    self.u_button.destroy()
                except:
                    pass
                
        for tab in self.tabs:
            if tab == self.tabs[0]:
                tabmenu.add_radiobutton(label=tab.title, command=(lambda tab=tab : set_value(tab)))
            else:
                tabmenu.add_radiobutton(label=tab.title, command=(lambda tab=tab : set_value(tab)))
        
        if len(self.tabs) > 1:
            self.menubar.add_cascade(label='Menus', menu=tabmenu)
        
        def show_help():
            tkMessageBox.showinfo(title='Help', message='To switch between' +\
            ' tabs, go to the Tabs menubar and select the desired tab')
        def show_about():
            tkMessageBox.showinfo(title='About', message=self.credits + 
            '\n\nVKML GUI created by: Alex Bartol\nCopyright 2009\nPurdue University' + \
            '\nNanohub.org')
        helpmenu = tk.Menu(self.menubar, tearoff=0)
        helpmenu.add_command(label='Help', command=show_help)
        helpmenu.add_command(label='About', command=show_about)
        self.menubar.add_cascade(label='Help', menu=helpmenu)
        
        self.top.config(menu=self.menubar)
        
    def get_result(self):
        """*Waits for the GUI window to close* and sets the value of 
        <Parameter>.out (or <Parameter>() ) to the desired value
        
        All illogical values are set to default
        """
        self.window.wait_window(self.top) #wait on window to close
        for i in range(len(self.results)):
            if type(self.args[i].type) != list and self.args[i].type != file:
                if self.args[i].withinrange(self.args[i].type(self.results[i])):
                    self.results[i] = self.args[i].type(self.results[i])
            else:
                self.results[i] = self.args[i].out
        i = 0
        for arg in self.args:
            arg.out = self.results[i]
            i += 1
    
    def _setup_gui(self):
        """sets up the GUI based on the number and name of each
        parameters
        
        NOTE: All variable types are handled uniquely
        """
        
        #tk.Label(self.top, text='Variable').grid(row=1, column=0, sticky=tk.W)
        #tk.Label(self.top, text='Value').grid(row=1, column=1, sticky=tk.W)

        self.top.bind("<Control-q>", self.close)
        #self.top.bind("<Button-1>", self.close)
        self.boolstatus.append(False) #needed to adjust for label row
        i = self.tabs[0].show(self.top)
        
        # if i == -1:
        #     pass
        # else:
        #     self.u_button = tk.Button(self.top, text=self.update_text,
        #                                 command=self._update_button)
        #     self.u_button.grid(row=i+3, sticky=tk.S,
        #                        column=0, columnspan=3)
        #     self.u_button.bind('<Return>', self._update_button)
        
        # c_button = tk.Button(self.top, text='Close',
        #    command=self.close)
        # c_button.grid(row=i+3, 
        #        column = 2, columnspan=3)
        # c_button.bind('<Return>', self.close)
    
    def _update(self):
        """This code is run whenever anything changes in the GUI
        
        It forwards the updates to the parser's callback function which runs
        each respective callback function if necessary
        """
        i = 0
        if len(self.results) == 0:
            for arg in self.args:
                self.results.append( arg.type(arg.get_widget_value()) )
        else:
            for i in range(len(self.args)):
                self.results[i] = self.args[i].type(
                                    self.args[i].get_widget_value())
        
        self.parser.update() #runs callbacks
        
    def _update_button(self, e=None):
        """This is used to bypass the interactive mode and force an update
        
        The callback functions are also informed of the update
        """
        i = 0
        if len(self.results) == 0:
            for arg in self.args:
                self.results.append( arg.type(arg.get_widget_value()) )
					  
        else:
            for i in range(len(self.args)):
                self.results[i] = self.args[i].type(
                                    self.args[i].get_widget_value())
				#print self.results[i]                    
        for i,x in enumerate(self.results):
            self.args[i].out = x
        self.parser.update_button() #forces an update on callback functions
        
    def close(self, e=None):
        """This method is run when the GUI should be destroyed
        """
        if self.parser.ask_quit:
            if tkMessageBox.askyesno(title='Quit Confirmation', 
                                    message='Do you really want to quit?') == 1:
                self.parser.run_post_commands()
                self.top.destroy()
                self.window.destroy()
                return
            return
        else:
            self.parser.run_post_commands()
            self.top.destroy()
            self.window.destroy()


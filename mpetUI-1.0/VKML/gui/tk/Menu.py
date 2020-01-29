"""Menu.py
"""
#
############################################################################
# FILE: Menu.py


#  Author: Alex Bartol <nanohub@alexbartol.com>
#       Copyright 2009
#
# 
try:
  #python 3
  import tkinter as tk
except ImportError:
  #python 2
  import Tkinter as tk
  
class Menu:
    """Menu is a class that holds parameters
    
    The Menu class should be used to organize <Parameter>s in any way the user
    see's fit
    """
    def __init__(self, parser, title='New Tab', update_button=True):
        """The <Menu> needs the parser because it automatically adds itself to
        it.
        
        parser - the parser the program is using
        title - the desired title for the Menu
        """
        self.update_button=update_button
        self.Parser = parser
        self.title = title
        self.params = []
        self.entries = []
        self.labels = []
        self.Parser.add_menu(self)

    def add(self, param):
        """This function adds param to the Menu
        
        param - a <Parameter> type
        """
        self.params.append(param)
        self.Parser.params.append(param)
        
    def show(self, top):
        """This function shows this Menu
        
        The Menu can be either shown or hidden depending on the user's current
        selection
        """
        for i, param in enumerate(self.params):
            self.labels.append(tk.Label(top, text=param.display_name))
            self.labels[-1].grid(row=i+3, column=0)
            self.entries.append(param.get_widget(top, row_number=i+3))
        if self.update_button:
	        return i+1
        return -1
        
    def hide(self):
        """This function hides the menu from the GUI
        
        Hiding a menu effectively temporarily destorys each parameter widget
        """
        for e in self.entries:
            e.destroy()
        for l in self.labels:
            l.destroy()
        for p in self.params:
            p.kill_special()
        self.entries = []
        self.labels = []
    

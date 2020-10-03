import tkinter as tk
from tkinter import *
from tkinter import ttk


class GUI(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        '''
        -------------------------------
        Frames for the UI:
            1. Labels Frame
            2. Goals Frame
            3. Molecule Frame
            4. Output Frame
        -------------------------------
        '''
        self.label_frame = tk.LabelFrame(self)
        self.label_frame.grid(row=0)

        self.goal_frame = tk.LabelFrame(self)
        self.goal_frame.grid(row=1)

        self.molecule_frame = tk.LabelFrame(self, width=450, height=200)
        self.molecule_frame.grid(row=2)

        self.output_frame = tk.LabelFrame(self)
        self.output_frame.grid(row=3)

        heading = tk.Label(self.label_frame, text="Training Environment for Drug Design", fg='gold')
        heading.grid(row=0)

        goalDropDownLabel = tk.Label(self.goal_frame, text="Optimisation goals")
        goalDropDownLabel.grid(row=0, column=1)

        # Iterations Label
        global counter
        counter = 0
        counterlabel = tk.Label(self.output_frame, text="Agent Iterations")
        counterlabel.grid(row=0)

        # Iterations Counter
        counter = tk.Label(self.output_frame, text=counter)
        counter.grid(row=1)

        # Preset Goal
        goalvalue = tk.StringVar()
        goal_chosen = ttk.Combobox(self.goal_frame, width=20, textvariable=goalvalue)
        goal_chosen["values"] = ("H20", "NaCl", "Minecraft Lava", "C3P0")
        goal_chosen.grid(column=1, row=1)
        goal_chosen.current(0)

        # Goal Button
        self.goalButton = tk.Button(self.goal_frame, text="START")
        self.goalButton["command"] = self.setgoal
        self.goalButton.grid(column=2, row=1)
        # self.goalButton.pack()

        # Quit Button
        quit = tk.Button(self.output_frame, text="QUIT", bg='black', command=self.master.destroy)
        quit.grid(row=2)
        # self.quit.grid(side="bottom")

        # Restart Button
        restart = tk.Button(self.goal_frame, text="RESTART", state='disabled')
        restart.grid(column=3, row=1)

    def setgoal(self):
        self.goalButton["text"] = "STOP"


root = tk.Tk()
root.geometry('1000x800')
# root.resizable(False, False)

root.title("TEDD.exe")
app = GUI(master=root)
app.mainloop()

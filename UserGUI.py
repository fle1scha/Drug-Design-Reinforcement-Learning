import tkinter as tk
from tkinter import ttk
from tkinter import messagebox as msg


class GUI(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_variables()
        self.create_widgets()

        self.paused = True

    def create_variables(self):
        self.counter = tk.IntVar()
        self.agentReward = tk.IntVar()
        self.agentAction = tk.IntVar()
        self.df = ("H20", "CO2", "NaCl", "C3PO")

    def create_widgets(self):
        """
        -------------------------------
        Frames for the UI:
            1. Labels Frame
            2. Goals Frame
            3. Molecule Frame
            4. Output Frame
        -------------------------------
        """

        # 1. Label Frame
        # ----------------------------
        self.label_frame = tk.LabelFrame(self)
        self.label_frame.grid(row=0)

        heading = tk.Label(self.label_frame, text="Training Environment for Drug Design")
        heading.grid(row=0)

        # ----------------------------

        # 2. Goals Frame
        # ----------------------------

        self.goal_frame = tk.LabelFrame(self, bg = "blue")
        self.goal_frame.grid(row=1)

        goalDropDownLabel = tk.Label(self.goal_frame, text="Optimisation goals")
        goalDropDownLabel.grid(row=0, column=1)

        # Preset Goal
        self.goalvalue = tk.StringVar()
        self.goal_chosen = ttk.Combobox(self.goal_frame, width=20, textvariable=self.goalvalue)
        self.goal_chosen["values"] = self.df
        self.goal_chosen.grid(column=1, row=1)
        self.goal_chosen.current(0)

        # Goal Button
        self.goalButton = tk.Button(self.goal_frame, text="START")
        self.goalButton["command"] = self.setgoal
        self.goalButton.grid(column=2, row=1)
        # self.goalButton.pack()

        # Restart Button
        self.restartButton = tk.Button(self.goal_frame, text="RESET", state='disabled', command=self.restart)
        self.restartButton.grid(column=3, row=1)

        # --------------------------

        # 3. Molecule Frame
        # -------------------------
        self.molecule_frame = tk.LabelFrame(self, width=450, height=200)
        self.molecule_frame.grid(row=2)

        self.output_frame = tk.LabelFrame(self)
        self.output_frame.grid(row=3)

        # ------------------------

        # 4. Output Frame
        # ------------------------

        # Iterations Label
        counterlabel = tk.Label(self.output_frame, text="Agent Iterations")
        counterlabel.grid(row=0, column=0)

        # Iterations Counter
        counter = tk.Label(self.output_frame, textvariable=self.counter)
        counter.grid(row=1, column=0)

        actionLabel = tk.Label(self.output_frame, text="Action Selected")
        actionLabel.grid(row=0, column=1)

        self.action = tk.Label(self.output_frame, textvariable=self.agentAction)
        self.action.grid(row=1, column=1)

        rewardLabel = tk.Label(self.output_frame, text="Recent Reward")
        rewardLabel.grid(row=0, column=2)

        reward = tk.Label(self.output_frame, textvariable=self.agentReward)
        reward.grid(row=1, column=2)

        # Quit Button
        quit = tk.Button(self.output_frame, text="QUIT", bg='black', command=self.master.destroy)
        quit.grid(row=3, column=1)

        # ------------------------

    """
    -------------------
    Button commands
    -------------------
    """

    def setgoal(self):

        if self.goalButton["text"] == "START":



            if self.goalvalue.get() == "":

                msg.showinfo("TEDD State Information", "You must enter an optimisation goal.")
            else:
                self.goal_chosen["state"] = 'disabled'
                self.restartButton["state"] = 'normal'
                self.goalButton["text"] = "PAUSE"
                return self.goalvalue.get()

        elif self.goalButton["text"] == "PAUSE":
            self.paused = True
            self.goalButton["text"] = "RESUME"

        elif self.goalButton["text"] == "RESUME":
            self.paused = False
            self.goalButton["text"] = "PAUSE"

    def restart(self):
        self.paused = True
        answer = msg.askyesno('TEDD State Information', "Are you sure you wish to reset the environment?")

        if answer:
            self.goal_chosen["state"] = 'normal'
            self.restartButton['state'] = 'disabled'
            answer = msg.askyesno('TEDD State Information', "Do you want to save the environment's lifespan to your "
                                                            "directory?")

            if answer:
                msg.showinfo('TEDD State Information', "The environment's lifespan has been saved in your directory.")
                self.goalButton["text"] = "START"
                self.goalvalue = ""
                return True
            if not answer:
                msg.showinfo('TEDD State Information', "The environment has reset.")
                self.goalButton["text"] = "START"
                self.goalvalue = ""

                return False

        if not answer:
            self.paused = False

        # publish list

    def invalidMol(self):
        msg.showerror('TEDD State Information', "You have entered an invalid molecule.")

        self.goalvalue = ""
        self.goalButton["text"] = "START"
        self.restartButton["state"] = "disabled"

    def get_counter(self):
        return self.counter

    def set_counter(self, value):
        self.counter.set(value)

    def get_agentreward(self):
        return self.agentReward

    def set_agentreward(self, value):
        self.agentReward.set(value)

    def get_agentaction(self):
        return self.agentAction

    def set_agentaction(self, value):
        self.agentAction.set(value)


    def getPaused(self):
        return self.state

    def setLibrary(self, df):
        self.df = df



#Creat a new GUI
root = tk.Tk()
root.geometry('500x400')
root.resizable(False, False)

root.title("TEDD.exe")

app = GUI(master=root)
app.goal_chosen.focus()
app.mainloop()

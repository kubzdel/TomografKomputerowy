import numpy
import tkinter
from PIL import Image, ImageTk
import tkinter as tk
from tkinter import Scale, Frame, Canvas
from tkinter import DoubleVar
from scipy import misc
import numpy as np
from random import random
data = np.zeros((200, 200))


class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        self.create_window()
        self.create_widgets()

    def create_widgets(self):
        self.hi_there = tk.Button(self)
        self.hi_there["text"] = "Hello World\n(click me)"
        self.hi_there["command"] = self.say_hi
        self.hi_there.pack(side="top")

        self.quit = tk.Button(self, text="QUIT", fg="red",
                              command=root.destroy)
        self.quit.pack(side="bottom")


    def create_window(self):
        self.frame = Frame(self.master, width=800, height=800)
        self.frame.pack()

        #data = misc.imread("/Users/martachecinska/Desktop/Kuba/semestrVI/InformatykaWMedycynie/Laboratoria/TomografKomputerowy/res/photo.png", flatten=True).astype('float64')
        self.canvas = Canvas(self.frame, width=len(data),height=len(data))
        self.canvas.place(x=-2,y=-2)
        self.im=Image.frombytes('L', (data.shape[1],data.shape[0]), data.astype('b').tostring())
        self.photo = ImageTk.PhotoImage(image=self.im)
        self.canvas.create_image(0,0,image=self.photo,anchor=tkinter.NW)
        self.master.update()

    def say_hi(self):
        print("hi there, everyone!")
    def random_image(self, value):
        data[(int)(random() * value), (int)(random()*value)] = 255
        self.im=Image.frombytes('L', (data.shape[1],data.shape[0]), data.astype('b').tostring())
        self.photo = ImageTk.PhotoImage(image=self.im)
        self.canvas.create_image(0,0,image=self.photo,anchor=tkinter.NW)
        self.canvas.update()
        self.master.update()

def method(value):
    app.random_image((int)(value))
root = tk.Tk()

var=DoubleVar()
app = Application(master=root)
w1 = Scale(root, from_=0, to_=150, tickinterval = 1, variable = var, command = method)
w1.set(45)
w1.pack()
w1.command = method
logo = tk.PhotoImage()

app.mainloop()
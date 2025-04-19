import tkinter as tk
from tkinter import messagebox


#Basic window to take user input
def main():
    print("")
    #Add functionality

#Set up window
window = tk.Tk()
window.title("Best Bus Route")
window.geometry("300x215")

#Starting location input box
startLabel = tk.Label(window, text="Enter starting location:")
startLabel.grid(row=0, column=0, padx=5, pady=5)
startBox = tk.Entry(window)
startBox.grid(row=0, column=1, padx=5, pady=5)

#Destination input box
endLabel = tk.Label(window, text="Enter destination:")
endLabel.grid(row=1, column=0, padx=5, pady=5)
endBox = tk.Entry(window)
endBox.grid(row=1, column=1, padx=5, pady=5)

#Use Dijkstra's toggle
toggleDij = tk.BooleanVar()
toggleDij_checkbox = tk.Checkbutton(window, text="Use Dijkstra's", variable=toggleDij)
toggleDij_checkbox.grid(row=2, columnspan=2, padx=5, pady=5)

#Use A* toggle
toggleA = tk.BooleanVar()
toggleA_checkbox = tk.Checkbutton(window, text="Use A*", variable=toggleA)
toggleA_checkbox.grid(row=3, columnspan=2, padx=5, pady=5)

#Toggle testing mode
test = tk.BooleanVar()
testBox = tk.Checkbutton(window, text="Testing Mode", variable=test)
testBox.grid(row=4, columnspan=2, padx=5, pady=5)

button = tk.Button(window, text="Find Route", command=main)
button.grid(row=5, columnspan=2, padx=5, pady=5)

window.mainloop()

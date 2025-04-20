import tkinter as tk
from tkinter import messagebox
from pathfinder import get_final_graph, getNodeFromLongLat, astar_shortest_path, dijkstra_shortest_path, visualize_path
import time
import webbrowser
def main():
    start = startBox.get()
    end = endBox.get()
    use_dij = toggleDij.get()
    use_astar = toggleA.get()
    testing = test.get()

    if not start or not end:
        messagebox.showerror("Missing Input", "Please enter both a starting location and a destination.")
        return

    if not use_dij and not use_astar:
        messagebox.showerror("No Algorithm Selected", "Please select at least one algorithm to use.")
        return

    try:
        G = get_final_graph()
        start_node, end_node = getNodeFromLongLat(start, end, G)

        results = []
        a_path = d_path = None

        # A*
        if use_astar:
            t0 = time.time()
            a_path, a_length = astar_shortest_path(G, start_node, end_node)
            a_time = time.time() - t0
            results.append(f"A* Path: {a_length:.2f} meters\nTime: {a_time:.4f} sec")

        # Dijkstra
        if use_dij:
            t0 = time.time()
            d_path, d_length = dijkstra_shortest_path(G, start_node, end_node)
            d_time = time.time() - t0
            results.append(f"Dijkstra Path: {d_length:.2f} meters\nTime: {d_time:.4f} sec")

        # Visualize and open
        filename = visualize_path(G, a_path if use_astar else None, d_path if use_dij else None)
        summary = "\n\n".join(results)
        summary += f"\n\nMap saved to '{filename}'."
        messagebox.showinfo("Route Found", summary)
        webbrowser.open(filename)

    except Exception as e:
        messagebox.showerror("Error", f"Something went wrong:\n{str(e)}")

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

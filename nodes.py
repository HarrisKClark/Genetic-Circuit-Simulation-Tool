import tkinter as tk

class Node:
    def __init__(self, canvas, x, y, node_id, radius=20, node_type='gene'):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.id = node_id
        self.radius = radius
        self.type = node_type

        fill_color = 'lightblue' if self.type == 'inducer' else 'white'

        self.oval = canvas.create_oval(
            x - radius, y - radius, x + radius, y + radius,
            fill=fill_color, outline='black', width=2,
            tags=('node', f'node{node_id}')
        )
        self.label = canvas.create_text(
            x, y, text=str(node_id),
            tags=('node_label', f'node{node_id}')
        )

    def set_color(self, color):
        self.canvas.itemconfig(self.oval, fill=color)

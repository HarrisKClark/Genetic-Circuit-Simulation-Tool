import tkinter as tk
import math

class Edge:
    def __init__(self, canvas, source, target, edge_type):
        self.canvas = canvas
        self.source = source
        self.target = target
        self.type = edge_type
        self.draw()

    def draw(self):
        if self.source is self.target:
            x, y = self.source.x, self.source.y
            r = self.source.radius
            loop_out = r * 2.0

            theta1 = math.radians(-45)
            theta2 = math.radians(45)

            start = (x + math.cos(theta1) * r, y - math.sin(theta1) * r)
            p1 = (x + math.cos(theta1) * (r + loop_out), y - math.sin(theta1) * (r + loop_out))
            p2 = (x + math.cos(theta2) * (r + loop_out), y - math.sin(theta2) * (r + loop_out))
            end = (x + math.cos(theta2) * r, y - math.sin(theta2) * r)

            if self.type == 'activator':
                self.canvas.create_line(
                    start[0], start[1],
                    p1[0], p1[1],
                    p2[0], p2[1],
                    end[0], end[1],
                    smooth=True, splinesteps=36,
                    width=2,
                    arrow=tk.LAST,
                    arrowshape=(12, 16, 6),
                    tags=('edge',)
                )
            else:
                self.canvas.create_line(
                    start[0], start[1],
                    p1[0], p1[1],
                    p2[0], p2[1],
                    end[0], end[1],
                    smooth=True, splinesteps=36,
                    width=2,
                    tags=('edge',)
                )
                tx = end[0] - p2[0]
                ty = end[1] - p2[1]
                tlen = math.hypot(tx, ty) or 1
                tx, ty = tx / tlen, ty / tlen
                px, py = -ty, tx
                bar = r * 0.7
                pA = (end[0] + px * bar / 2, end[1] + py * bar / 2)
                pB = (end[0] - px * bar / 2, end[1] - py * bar / 2)
                self.canvas.create_line(
                    pA[0], pA[1],
                    pB[0], pB[1],
                    width=2,
                    tags=('edge',)
                )
            return

        x1, y1 = self.source.x, self.source.y
        x2, y2 = self.target.x, self.target.y
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy)
        if dist == 0:
            return

        ux, uy = dx / dist, dy / dist
        start = (x1 + ux * self.source.radius, y1 + uy * self.source.radius)
        end = (x2 - ux * self.target.radius, y2 - uy * self.target.radius)

        mx, my = (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
        perp_x, perp_y = -uy, ux
        curve_off = min(100, dist / 2)
        ctrl = (mx + perp_x * curve_off, my + perp_y * curve_off)

        if self.type == 'activator':
            self.canvas.create_line(
                start[0], start[1],
                ctrl[0], ctrl[1],
                end[0], end[1],
                smooth=True, splinesteps=20,
                width=2,
                arrow=tk.LAST,
                arrowshape=(12, 16, 6),
                tags=('edge',)
            )
        else:
            self.canvas.create_line(
                start[0], start[1],
                ctrl[0], ctrl[1],
                end[0], end[1],
                smooth=True, splinesteps=20,
                width=2,
                tags=('edge',)
            )
            tx, ty = end[0] - ctrl[0], end[1] - ctrl[1]
            tlen = math.hypot(tx, ty) or 1
            tx, ty = tx / tlen, ty / tlen
            px, py = -ty, tx
            bar = 16
            p1 = (end[0] + px * bar / 2, end[1] + py * bar / 2)
            p2 = (end[0] - px * bar / 2, end[1] - py * bar / 2)
            self.canvas.create_line(
                p1[0], p1[1],
                p2[0], p2[1],
                width=2,
                tags=('edge',)
            )

import numpy as np
from scipy.integrate import solve_ivp
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.animation as animation

class ThreeBodyProblem:
    def __init__(self, masses, initial_conditions, time_span):
        self.G = 1.0  # Gravitational constant
        self.masses = masses
        self.initial_conditions = initial_conditions
        self.time_span = time_span

    def equations_of_motion(self, t, state):
        x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 = state
        
        r12 = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        r13 = np.sqrt((x3 - x1)**2 + (y3 - y1)**2)
        r23 = np.sqrt((x3 - x2)**2 + (y3 - y2)**2)
        
        ax1 = self.G * (self.masses[1] * (x2 - x1) / r12**3 + self.masses[2] * (x3 - x1) / r13**3)
        ay1 = self.G * (self.masses[1] * (y2 - y1) / r12**3 + self.masses[2] * (y3 - y1) / r13**3)
        ax2 = self.G * (self.masses[0] * (x1 - x2) / r12**3 + self.masses[2] * (x3 - x2) / r23**3)
        ay2 = self.G * (self.masses[0] * (y1 - y2) / r12**3 + self.masses[2] * (y3 - y2) / r23**3)
        ax3 = self.G * (self.masses[0] * (x1 - x3) / r13**3 + self.masses[1] * (x2 - x3) / r23**3)
        ay3 = self.G * (self.masses[0] * (y1 - y3) / r13**3 + self.masses[1] * (y2 - y3) / r23**3)

        return [vx1, vy1, ax1, ay1, vx2, vy2, ax2, ay2, vx3, vy3, ax3, ay3]

    def solve(self, time_points):
        solution = solve_ivp(self.equations_of_motion, self.time_span, self.initial_conditions, 
                             t_eval=time_points, method='RK45', rtol=1e-8, atol=1e-8)
        return solution

class ThreeBodyAnimator:
    def __init__(self, solution):
        self.solution = solution
        self.root = tk.Tk()
        self.root.title("Three-Body Problem Animation")

        self.fig = Figure(figsize=(10, 10), dpi=100)
        self.ax = self.fig.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.time_text = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)

    def init_animation(self):
        self.ax.clear()
        self.ax.set_xlim(min(self.solution.y[0].min(), self.solution.y[4].min(), self.solution.y[8].min()) - 0.5,
                         max(self.solution.y[0].max(), self.solution.y[4].max(), self.solution.y[8].max()) + 0.5)
        self.ax.set_ylim(min(self.solution.y[1].min(), self.solution.y[5].min(), self.solution.y[9].min()) - 0.5,
                         max(self.solution.y[1].max(), self.solution.y[5].max(), self.solution.y[9].max()) + 0.5)
        self.ax.set_xlabel('X Position')
        self.ax.set_ylabel('Y Position')
        self.ax.set_title('Three-Body Problem Animation')
        return []

    def update_animation(self, frame):
        self.ax.clear()
        self.init_animation()
        
        # Plot trajectories
        self.ax.plot(self.solution.y[0][:frame], self.solution.y[1][:frame], 'b-', lw=1, alpha=0.5)
        self.ax.plot(self.solution.y[4][:frame], self.solution.y[5][:frame], 'g-', lw=1, alpha=0.5)
        self.ax.plot(self.solution.y[8][:frame], self.solution.y[9][:frame], 'r-', lw=1, alpha=0.5)
        
        # Plot current positions
        self.ax.plot(self.solution.y[0][frame], self.solution.y[1][frame], 'bo', markersize=10, label='Body 1')
        self.ax.plot(self.solution.y[4][frame], self.solution.y[5][frame], 'go', markersize=10, label='Body 2')
        self.ax.plot(self.solution.y[8][frame], self.solution.y[9][frame], 'ro', markersize=10, label='Body 3')
        
        self.ax.legend()
        self.time_text.set_text(f'Time: {self.solution.t[frame]:.2f}')
        self.canvas.draw()

    def animate(self):
        anim = animation.FuncAnimation(self.fig, self.update_animation, frames=len(self.solution.t),
                                       interval=50, repeat=False)
        self.root.mainloop()

def main():
    # Simulation parameters
    masses = [1.0, 1.0, 1.0]
    initial_conditions = [0, 1, 0.5, 0, 1, 0, 0, 0.5, -1, -1, -0.5, 0]
    time_span = (0, 50)
    time_points = np.linspace(0, 50, 500)

    # Solve the three-body problem
    print("Solving the three-body problem...")
    three_body = ThreeBodyProblem(masses, initial_conditions, time_span)
    solution = three_body.solve(time_points)

    # Create and display the animation
    print("Generating and displaying animation...")
    animator = ThreeBodyAnimator(solution)
    animator.animate()

if __name__ == "__main__":
    main()
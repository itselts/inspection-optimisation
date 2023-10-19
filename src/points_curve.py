import numpy as np
import matplotlib.pyplot as plt

def stage_one(x, growth, due_score, minimum_score):
    '''Period before due date.'''
    return (due_score - minimum_score) * growth ** (x) + minimum_score

def stage_two(x, height, maxi, width):
    '''Single cycle immediately after due date.'''
    return (height + 1) ** (((x % width))/width) + maxi - (height + 1)

def stage_three(x, height, maxi, width):
    '''Repeated cycle after due date.'''
    return (height + 1) ** (((x % width))/width) + maxi - (height + 1)

def combined(x, growth, due_score, minimum_score, cycle_height1, cycle_maximum1, width1, cycle_height2, cycle_maximum2, width2, act_multi):
    '''Combining all three stages.'''
    result = np.piecewise(
        x, 
        [x <= 0, (x>0)*(x<10), x >= 10], 
        [
            lambda x : stage_one(x, growth, due_score, minimum_score), 
            lambda x : stage_two(x, cycle_height1, cycle_maximum1, width1),
            lambda x : stage_three(x, cycle_height2, cycle_maximum2, width2)
        ]
    )
    return result * act_multi

def plot_curve(request_ID, x = np.linspace(-30, 30, 10000)):
    '''Plots the points curve and saves it to the results folder.'''
    fig1 = plt.figure(0)
    ax = fig1.add_subplot(1, 1, 1)
    ax.grid()
    line, = ax.plot(x, combined(
        x = x, 
        growth = 1.25,
        due_score = 250,
        minimum_score = 50,
        cycle_height1 = 150, 
        cycle_maximum1 = 250,
        width1 = 10,
        cycle_height2 = 100, 
        cycle_maximum2 = 250,
        width2 = 5,
        act_multi = 1))
    ax.set_title("Points curve based on due date")
    ax.set_xlabel("Negative days until due date")
    ax.set_ylabel("Points")

    # Plotting the points curve and saving it
    fig1.savefig(f'../runs/{request_ID}/outputs/points_curve.png')
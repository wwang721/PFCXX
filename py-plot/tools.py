# Some math tools.
# Created by W. Wang

import numpy as np
import matplotlib.pyplot as plt
import imageio


def rfKMC_random():
    """ Uniform random number in (0,1].
    """
    u = np.random.rand()
    if u == 0.:
        u = 1.
    return u


def binary_search(target, source):
    """ Binary search the position of target in the source set.
    """
    N = len(source)
    start = 0
    end = N-1

    while (end - start) != 1:
        mid = start + (end - start)//2

        if target <= source[mid]:
            end = mid
        else:
            start = mid

    return start


def save_ani(image_list, file_name, duration=0.35, fps=5):
    """ Saves the animation created.
    """
    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))

    #imageio.mimsave(file_name+".gif", frames, 'GIF', duration=duration)
    imageio.mimsave(file_name+".mp4", frames, 'MP4', fps=fps)


def integration(x, y, yerr):
    """ Trapezoidal method to compute integral and corresponding error
    """
    N = len(x)
    Ans = 0.
    err = 0.
    for k in range(N-1):
        dx = (x[k+1]-x[k])
        Ans += (y[k]+y[k+1]) * dx/2
        deltaAns = np.sqrt(yerr[k]**2+yerr[k+1]**2)*dx/2
        err += deltaAns**2
    err = np.sqrt(err)
    return Ans, err


def main():
    # Test binary search function
    test_search = [0, 0, 0, 0, 0, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 4]
    print(binary_search(4, test_search))

    # Test single junction rfKMC
    k = 0.2  # failure rate
    T = []  # container for rupture time
    N = 1000  # number of assays
    for i in range(N):
        T.append(np.log(1/rfKMC_random())/k)

    T = np.array(T)  # transform to numpy array

    S = []  # container for survival curve
    tt = np.linspace(0, T.max(), 100)
    for t in tt:
        S.append(len(T[T > t])/N)

    fig, ax = plt.subplots()
    ax.plot(tt, S, label="rfKMC")
    P = [np.exp(-k*t) for t in tt]  # prediction of the survival curve
    ax.plot(tt, P, "--", label="Prediction")
    ax.set_xlim(0, tt[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"Time")
    ax.set_ylabel(r"Survival probability")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()

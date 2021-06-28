import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from tools import save_ani


def angle_regulate(theta):
    """ regulates angles > -pi and <= pi
    """
    theta = theta % (2*np.pi)
    while theta <= -np.pi:
        theta += 2 * np.pi
    while theta > np.pi:
        theta -= 2 * np.pi
    return theta


def chi(theta, x, y):
    """ Preferrence of angle
    """
    k = 0.01
    angle = angle_regulate(theta)

    if y > -10:
        preferred_angle = np.pi/2
    else:
        preferred_angle = np.arctan2(-10-y, -x)

    return k*(angle - preferred_angle)


def random_angle(x, y):
    """ Generates Gaussian random angle
    """
    D = 0.001
    k = 0.01
    sigma = np.sqrt(D/k)
    if y > -10:
        preferred_angle = np.pi/2
    else:
        preferred_angle = np.arctan2(-10-y, -x)
    return angle_regulate(np.random.normal(preferred_angle, sigma))


class field(object):

    def __init__(self, x0, y0, r0, v0=[0.025, 0.025]):

        self.position = (x0, y0)
        thickness = 2
        self.r0 = r0
        self.pixel = 1.
        self.v0 = v0
        self.vv = np.linalg.norm(v0)
        self.polarity = np.arctan2(v0[1], v0[0])

        # subdomain center
        self.subdomain = [np.floor(x0/self.pixel), np.floor(y0/self.pixel)]
        # center of cell in subdomain
        self.center = (
            x0 - self.subdomain[0] * self.pixel, y0 - self.subdomain[1] * self.pixel)

        self.N = 50
        self.phi = np.zeros((self.N+1, self.N+1))
        self.indices = [(i, j) for i in range(self.N+1)
                        for j in range(self.N+1)]
        self.phi0 = 1.

        for (i, j) in self.indices:
            pos = self.index2position(i, j)
            rvec = (pos[0]-self.center[0], pos[1]-self.center[1])
            r = np.linalg.norm(rvec)
            self.phi[i, j] = self.phi0/2 * \
                np.tanh((r0 - r)/thickness) + self.phi0/2

        self.fullN = (160, 100)
        self.fulllattice()

    def index2position(self, i, j):
        x = (j - self.N/2) * self.pixel
        y = (i - self.N/2) * self.pixel
        return (x, y)

    def subdomain_shift(self, x, y):
        jshift = int(np.floor(x/self.pixel))
        ishift = int(np.floor(y/self.pixel))
        temp = np.zeros(self.phi.shape)

        for (i, j) in self.indices:
            ii = i + ishift
            jj = j + jshift

            if ii > self.N:
                ii = ii - self.N - 1
            elif ii < 0:
                ii = self.N + ii + 1

            if jj > self.N:
                jj = jj - self.N - 1
            elif jj < 0:
                jj = self.N + jj + 1

            temp[i, j] = self.phi[ii, jj]
        self.phi = temp

        self.subdomain[0] += jshift
        self.subdomain[1] += ishift

        if self.subdomain[0] > self.fullN[1]/2:
            self.subdomain[0] = self.subdomain[0] - self.fullN[1]
        elif subdomain[0] < -fullN[1]/2:
            self.subdomain[0] = self.subdomain[0] + self.fullN[1]

        if self.subdomain[1] > self.fullN[0]/2:
            self.subdomain[1] = self.subdomain[1] - self.fullN[0]
        elif self.subdomain[1] < -self.fullN[0]/2:
            self.subdomain[1] = self.subdomain[1] + self.fullN[0]

    def shift(self):
        self.center = self.center_of_mass()
        d = np.linalg.norm(self.center)
        if d > 2 * self.pixel:
            self.subdomain_shift(self.center[0], self.center[1])

    def fulllattice(self):
        fullN = self.fullN
        self.fullphi = np.zeros((fullN[0]+1, fullN[1]+1))
        for (i, j) in self.indices:

            ii = int(i - self.N/2 + self.subdomain[1] + fullN[0]/2)
            jj = int(j - self.N/2 + self.subdomain[0] + fullN[1]/2)

            if ii > fullN[0]:
                ii = ii - fullN[0] - 1
            elif ii < 0:
                ii = fullN[0] + ii + 1

            if jj > fullN[1]:
                jj = jj - fullN[1] - 1
            elif jj < 0:
                jj = fullN[1] + jj + 1

            self.fullphi[ii, jj] = self.phi[i, j]

    def fulllatice2subdomain(self, h):
        fullN = self.fullN
        hs = np.zeros((self.N+1, self.N+1))
        for (i, j) in self.indices:

            ii = int(i - self.N/2 + self.subdomain[1] + fullN[0]/2)
            jj = int(j - self.N/2 + self.subdomain[0] + fullN[1]/2)

            if ii > fullN[0]:
                ii = ii - fullN[0] - 1
            elif ii < 0:
                ii = fullN[0] + ii + 1

            if jj > fullN[1]:
                jj = jj - fullN[1] - 1
            elif jj < 0:
                jj = fullN[1] + jj + 1

            hs[i, j] = h[ii, jj]

        return hs

    def central_diff(self, phi, axis):
        ''' periodic boundary condition
        '''
        # axis = 0: y-axis
        # axis = 1: x-axis
        (N1, N2) = phi.shape
        res = np.zeros((N1, N2))
        if axis == 0:
            for i in range(N1):
                for j in range(N2):
                    if i+1 >= N1:
                        res[i, j] = (phi[0, j] - phi[i-1, j])/(2*self.pixel)
                    elif i-1 < 0:
                        res[i, j] = (phi[i+1, j] - phi[N1-1, j])/(2*self.pixel)
                    else:
                        res[i, j] = (phi[i+1, j] - phi[i-1, j])/(2*self.pixel)
        elif axis == 1:
            for i in range(N1):
                for j in range(N2):
                    if j+1 >= N2:
                        res[i, j] = (phi[i, 0] - phi[i, j-1])/(2*self.pixel)
                    elif j-1 < 0:
                        res[i, j] = (phi[i, j+1] - phi[i, N2-1])/(2*self.pixel)
                    else:
                        res[i, j] = (phi[i, j+1] - phi[i, j-1])/(2*self.pixel)
        return res

    def second_diff(self, phi, axis):
        ''' periodic boundary condition
        '''
        # axis = 0: y-axis
        # axis = 1: x-axis

        (N1, N2) = phi.shape
        res = np.zeros((N1, N2))
        if axis == 0:
            for i in range(N1):
                for j in range(N2):
                    if i+1 >= N1:
                        res[i, j] = (phi[0, j] - 2*phi[i, j] +
                                     phi[i-1, j])/(self.pixel**2)
                    elif i-1 < 0:
                        res[i, j] = (phi[i+1, j] - 2*phi[i, j] +
                                     phi[N1-1, j])/(self.pixel**2)
                    else:
                        res[i, j] = (phi[i+1, j] - 2*phi[i, j] +
                                     phi[i-1, j])/(self.pixel**2)
        elif axis == 1:
            for i in range(N1):
                for j in range(N2):
                    if j+1 >= N2:
                        res[i, j] = (phi[i, 0] - 2*phi[i, j] +
                                     phi[i, j-1])/(self.pixel**2)
                    elif j-1 < 0:
                        res[i, j] = (phi[i, j+1] - 2*phi[i, j] +
                                     phi[i, N2-1])/(self.pixel**2)
                    else:
                        res[i, j] = (phi[i, j+1] - 2*phi[i, j] +
                                     phi[i, j-1])/(self.pixel**2)
        return res

    def area(self):
        A = 0.
        for i in range(self.N):
            for j in range(self.N):
                p = (self.phi[i, j] + self.phi[i+1, j] +
                     self.phi[i, j+1] + self.phi[i+1, j+1]) / 4
                A += (p**2) * (self.pixel**2)
        return A

    def center_of_mass(self):
        x = 0.
        y = 0.
        m = 0.
        for i in range(self.N):
            for j in range(self.N):
                (xx, yy) = self.index2position(i, j)
                p = (self.phi[i, j] + self.phi[i+1, j] +
                     self.phi[i, j+1] + self.phi[i+1, j+1]) / 4
                m += p * (self.pixel**2)
                x += p * xx * (self.pixel**2)
                y += p * yy * (self.pixel**2)
        return (x/m, y/m)

    def update_polarity(self, dt, chemo=False):
        mu = 0.  # mean
        D = 0.001  # diffusion coefficient
        sigma = np.sqrt(2*D*dt)  # standard deviation
        self.polarity += np.random.normal(mu, sigma)

        if chemo is True:
            x = self.center[0] + self.subdomain[0] * self.pixel
            y = self.center[1] + self.subdomain[1] * self.pixel
            self.polarity -= chi(self.polarity, x, y) * dt

        self.v0 = [self.vv * np.cos(self.polarity),
                   self.vv * np.sin(self.polarity)]

    def update(self, dt=0.5, h=np.zeros((10, 10)), h1_laplace=np.zeros((10, 10)), c=np.zeros((100, 100)),
               multi=False, adh=False, confine=False, chemo=False):

        g_x = self.central_diff(self.phi, 1)
        g_y = self.central_diff(self.phi, 0)
        laplace = self.second_diff(self.phi, 0) + self.second_diff(self.phi, 1)
        deltaV = 1. - (self.area() / (np.pi * ((self.r0 * self.phi0)**2)))

        if multi is True:
            hs = self.fulllatice2subdomain(h)

        if confine is True:
            cs = self.fulllatice2subdomain(c)
            epsilon_c = 1.

        if adh is True:
            h1_laplace_sub = self.fulllatice2subdomain(h1_laplace)
            omega = 0.05

        alpha = 1.
        K = 2.
        lam = 600
        gamma = 10
        epsilon = 0.2

        for (i, j) in self.indices:
            dphi1 = -np.dot(self.v0, [g_x[i, j], g_y[i, j]])
            dphi2 = alpha * (self.phi[i, j]**3 - 1.5 * (self.phi[i, j]**2)
                             * self.phi0 + 0.5 * self.phi[i, j] * (self.phi0**2))
            dphi3 = - K * laplace[i, j]
            dphi4 = - \
                (4 * lam * self.phi[i, j] / (np.pi *
                                             ((self.r0 * self.phi0)**2))) * deltaV
            if multi is True:
                dphi5 = 2 * epsilon * \
                    self.phi[i, j] * (hs[i, j] - (self.phi[i, j])**2)
            else:
                dphi5 = 0.

            if confine is True:
                dphi6 = 2 * epsilon_c * self.phi[i, j] * (cs[i, j]**2)
            else:
                dphi6 = 0.

            if adh is True:
                dphi_adh = - omega * (h1_laplace_sub[i, j] - laplace[i, j])
            else:
                dphi_adh = 0.

            self.phi[i, j] += (dphi1 - ((dphi2 + dphi3 +
                                         dphi4 + dphi5 + dphi6 + dphi_adh)/gamma)) * dt

        self.shift()
        self.fulllattice()
        self.update_polarity(dt, chemo)

    def showpolarity(self):
        scale = 250
        x0 = self.center[0] + self.subdomain[0] * self.pixel
        y0 = self.center[1] + self.subdomain[1] * self.pixel

        plt.arrow(x0, y0, self.v0[0] * scale, self.v0[1] * scale,
                  width=0.8,
                  length_includes_head=True,  # 增加的长度包含箭头部分
                  head_width=2.8,
                  head_length=4,
                  fc=(0/255, 52/255, 114/255),
                  ec=(0/255, 52/255, 114/255))

    def showdomainbox(self):
        x = (self.subdomain[0] - self.N/2) * self.pixel
        y = (self.subdomain[1] - self.N/2) * self.pixel
        plt.plot([x, x], [y, y+self.N * self.pixel])
        plt.plot([x+self.N * self.pixel, x + self.N * self.pixel],
                 [y, y+self.N * self.pixel])
        plt.plot([x, x+self.N * self.pixel], [y, y])
        plt.plot([x, x+self.N * self.pixel],
                 [y+self.N * self.pixel, y+self.N * self.pixel])

    def showplot(self, t):
        edge0 = self.fullN[0] * self.pixel
        edge1 = self.fullN[1] * self.pixel

        x = np.linspace(-edge1/2, edge1/2, self.fullN[1]+1)
        y = np.linspace(-edge0/2, edge0/2, self.fullN[0]+1)

        X, Y = np.meshgrid(x, y)

        self.showdomainbox()
        self.showpolarity()
        plt.imshow(self.fullphi, cmap='gray_r', origin="lower", extent=(-edge1/2, edge1/2, -
                                                                        edge0/2, edge0/2), vmax=1.5, vmin=0)		# plt.contour(X,Y, phi, [phi0/2], colors="k")
        plt.colorbar()
        plt.contour(X, Y, self.fullphi, [self.phi0/2], colors="k")
        plt.savefig("figures/%s.png" % t)
        # plt.show()

    def simulation(self, T, dt, ti=0.):
        t = ti
        k = 0
        plt.figure()
        self.showplot(t)
        # plt.cla()

        pbar = tqdm(total=T//dt, desc="Processing")
        while t <= T+ti:
            k += 1
            t += dt
            self.update(dt)
            if k % 10 == 0:
                plt.clf()
                self.showplot(t)
            # plt.pause(1)
            pbar.update(1)
            # print(t)

    def animation(self, T, dt, file_name):
        t = 0.
        k = 0
        image_list = ["figures/%s.png" % t]
        while t <= T:
            k += 1
            t += dt
            if k % 10 == 0:
                image_list.append("figures/%s.png" % t)

        duration = 0.2
        fps = 5
        save_ani(image_list, file_name, duration, fps)


def main():
    a = field(0, 0, 12)
    # a.showplot(0)

    a.simulation(400, 0.5)
    a.animation(400, 0.5, "test")
    '''
    edge = 50.
    plt.figure()
    plt.imshow(a.second_diff(0), cmap='gray_r', origin="lower", extent=(-edge/2, edge/2, -edge/2, edge/2))
    plt.colorbar()
    plt.show()
    '''


if __name__ == "__main__":
    main()

import numpy as np
import matplotlib.pyplot as plt
from single_field import field
from tqdm import tqdm
from tools import rfKMC_random
import h5py


class system(field):

    def __init__(self):
        self.sys = []
        x0s = [-40, -20, 0, 20, 40]
        y0s = [-60, -40, -20]
        for x0 in x0s:
            for y0 in y0s:
                angle = (rfKMC_random()-0.5) * 2. * np.pi
                # angle = random_angle(x0, y0)
                v0 = [0.035 * np.cos(angle), 0.035 * np.sin(angle)]
                self.sys.append(field(x0, y0, 12, v0))

        self.fullN = self.sys[0].fullN
        self.pixel = self.sys[0].pixel

        self.confinement = np.zeros((self.fullN[0]+1, self.fullN[1]+1))
        self.confinement[70:, :40] = np.ones(self.confinement[70:, :40].shape)
        self.confinement[70:, 60:] = np.ones(self.confinement[70:, 60:].shape)
        self.confinement[:10, :] = np.ones(self.confinement[:10, :].shape)

        # for thermalization
        self.conf = np.zeros((self.fullN[0]+1, self.fullN[1]+1))
        self.conf[70:, :] = np.ones(self.conf[70:, :].shape)
        self.conf[:10, :] = np.ones(self.conf[:10, :].shape)

    def thermalization_update(self, dt=0.5):
        fullN = self.fullN
        h = np.zeros((fullN[0]+1, fullN[1]+1))  # for perimeter (tension) term
        h1 = np.zeros((fullN[0]+1, fullN[1]+1)) 	# for adhesion term
        for f in self.sys:
            h += f.fullphi**2
            h1 += f.fullphi

        h1_laplace = self.second_diff(h1, 0) + self.second_diff(h1, 1)

        for f in self.sys:
            f.update(dt, h, h1_laplace, self.conf,
                     multi=True, adh=True, confine=True, chemo=False)

    def thermalization(self, T, dt):
        t = 0.
        k = 0
        plt.figure()
        self.showplot(t, therm=True)
        # plt.cla()

        pbar = tqdm(total=T//dt, desc="Thermalize")
        while t <= T:
            k += 1
            t += dt
            self.thermalization_update(dt)
            if k % 10 == 0:
                plt.clf()
                self.showplot(t, therm=True)
            # plt.pause(1)
            pbar.update(1)

    def update(self, dt=0.5):
        fullN = self.fullN
        h = np.zeros((fullN[0]+1, fullN[1]+1))  # for perimeter (tension) term
        h1 = np.zeros((fullN[0]+1, fullN[1]+1)) 	# for adhesion term
        for f in self.sys:
            h += f.fullphi**2
            h1 += f.fullphi

        h1_laplace = self.second_diff(h1, 0) + self.second_diff(h1, 1)

        for f in self.sys:
            f.update(dt, h, h1_laplace, self.confinement,
                     multi=True, adh=True, confine=True, chemo=True)

    def showplot(self, t, fileName, therm=False):
        edge0 = self.fullN[0] * self.pixel
        edge1 = self.fullN[1] * self.pixel

        phi = np.zeros((self.fullN[0]+1, self.fullN[1]+1))
        '''
        for f in self.sys:
            phi += f.fullphi
        '''
        if therm is True:
            phi += self.conf
        else:
            phi += self.confinement

        x = np.linspace(-edge1/2, edge1/2, self.fullN[1]+1)
        y = np.linspace(-edge0/2, edge0/2, self.fullN[0]+1)

        X, Y = np.meshgrid(x, y)

        plt.imshow(phi, cmap='gray_r', origin="lower", extent=(-edge1 /
                                                               2, edge1/2, -edge0/2, edge0/2), vmax=1.5, vmin=0)
        # plt.imshow(self.sys[0].phi, cmap='gray_r', origin="lower", extent=(-edge/2, edge/2, -edge/2, edge/2), vmax=1.5, vmin=0)
        # plt.imshow(self.sys[1].phi, cmap='gray_r', origin="lower", extent=(-edge/2, edge/2, -edge/2, edge/2), vmax=1.5, vmin=0)		# plt.contour(X,Y, phi, [phi0/2], colors="k")
        # plt.colorbar()
        # plt.contour(X,Y, phi, [self.sys[0].phi0/2], colors="k")
        # self.sys[0].showdomainbox()
        for f in self.sys:
            f.showpolarity()
            # field.showdomainbox()
            plt.contour(X, Y, f.fullphi, [f.phi0/2], colors="k")
            plt.contourf(X, Y, f.fullphi, [f.phi0/2, (f.phi0/2) + 15*f.phi0], colors=[(127/255,236/255,173/255)], alpha=0.3)

        if therm is True:
            plt.contour(X, Y, self.conf, [0.5], colors="k")
        else:
            plt.contour(X, Y, self.confinement, [0.5], colors="k")

        plt.xlim(-edge1/2, edge1/2)
        plt.ylim(-edge0/2, edge0/2)
        plt.title("T=%g" % t)
        # plt.show()
        plt.savefig("figures/%s.png" % fileName, dpi=300)

    def save_data(self, file_name="data.h5"):
        file = h5py.File(file_name, "w")
        index = 0
        for f in self.sys:
            data = file.create_dataset(
                "dset{}".format(index), data=f.phi)
            data.attrs["subdomain"] = f.subdomain
            data.attrs["center"] = f.center
            data.attrs["velocity"] = f.v0
            index += 1
        file.close()

    def load_data(self, file_name="data.h5"):
        file = h5py.File(file_name, "r")
        index = 0
        for f in self.sys:
            dataset = "dset{}".format(index)
            f.phi = file[dataset][:]
            f.subdomain = file[dataset].attrs["subdomain"]
            f.center = file[dataset].attrs["center"]
            f.v0 = file[dataset].attrs["velocity"]
            
            f.fulllattice()
            f.vv = np.linalg.norm(f.v0)
            f.polarity = np.arctan2(f.v0[1], f.v0[0])
            index += 1

        file.close()


def main():
    sys = system()
    sys.thermalization(10, 0.5)
    sys.save_data()
    sys.load_data()
    sys.simulation(10, 0.5, 10)
    sys.animation(20, 0.5, "test")


if __name__ == "__main__":
    main()

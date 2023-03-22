# -----------------------------------------------------------------------------
# Space Constellation Simulation and Analysis
# -----------------------------------------------------------------------------
#
# Authors: Subramanian Arumugam (Astrodynmaics Researcher)
#
# Description:
# This Python script simulates a constellation of satellites and space debris in
# Earth orbit. It propagates their orbits over time and calculates the miss
# distances between satellites and debris. The script also saves the positions
# of the satellites and miss distances as separate CSV files.
#
# -----------------------------------------------------------------------------


import numpy as np
import csv
import os
from astropy import units as u
from astropy.time import Time
from poliastro.twobody import Orbit
from poliastro.bodies import Earth
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from astropy.coordinates import CartesianRepresentation

class SpaceObject:
    """Base class for space objects, like satellites and debris."""

    def __init__(self, semi_major_axis, eccentricity, inclination, raan, arg_perigee, true_anomaly, epoch):
        self.orbit = Orbit.from_classical(
            Earth,
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            arg_perigee,
            true_anomaly,
            epoch=epoch
        )

    def propagate(self, time):
        """Propagate the space object's orbit."""
        return self.orbit.propagate(time)

class Satellite(SpaceObject):
    """Satellite class, inheriting from SpaceObject."""

    def __init__(self, altitude, inclination, raan, true_anomaly):
        semi_major_axis = (Earth.R + altitude).to(u.m)
        super().__init__(semi_major_axis, 0 * u.one, inclination, raan, 90 * u.deg, true_anomaly, Time("J2000"))

class Debris(SpaceObject):
    """Space debris class, inheriting from SpaceObject."""

    def __init__(self, altitude, inclination, raan, true_anomaly):
        semi_major_axis = (Earth.R + altitude).to(u.m)
        super().__init__(semi_major_axis, 0 * u.one, inclination, raan, 90 * u.deg, true_anomaly, Time("J2000"))

class Constellation:
    """Class to create and manage a constellation of space objects."""

    def __init__(self):
        self.space_objects = []

    def add_space_object(self, obj):
        """Add a space object to the constellation."""
        self.space_objects.append(obj)

    def create_sso_constellation(self, N, altitude, num_planes):
        """Create a sun-synchronous orbit constellation."""
        raans = np.linspace(0, 360, num_planes + 1)[:-1] * u.deg
        true_anomalies = np.linspace(0, 360, N // num_planes + 1)[:-1] * u.deg

        for raan in raans:
            for true_anomaly in true_anomalies:
                satellite = Satellite(altitude, 98 * u.deg, raan, true_anomaly)
                self.add_space_object(satellite)

    def create_space_debris(self, N_debris, altitude, inclination_range):
        """Create space debris objects."""
        inclinations = np.random.uniform(inclination_range[0], inclination_range[1], N_debris) * u.deg
        raans = np.random.uniform(0, 360, N_debris) * u.deg
        true_anomalies = np.random.uniform(0, 360, N_debris) * u.deg

        for inclination, raan, true_anomaly in zip(inclinations, raans, true_anomalies):
            debris = Debris(altitude, inclination, raan, true_anomaly)
            self.add_space_object(debris)

    def calculate_miss_distances(self):
        """Calculate miss distances between satellites and debris objects."""
        satellites = [obj for obj in self.space_objects if isinstance(obj, Satellite)]
        debris_objects = [obj for obj in self.space_objects if isinstance(obj, Debris)]
        miss_distances = []

        for sat in satellites:
            sat_distances = []
            for deb in debris_objects:
                sat_pos = sat.orbit.r.to(u.km).value
                deb_pos = deb.orbit.r.to(u.km).value
                distance = np.linalg.norm(sat_pos - deb_pos)
                sat_distances.append(distance)
            miss_distances.append(sat_distances)

        return miss_distances

    def save_positions_and_distances_to_csv(self, time_steps):
        """Save satellite positions and miss distances to separate CSV files."""
        satellites = [obj for obj in self.space_objects if isinstance(obj, Satellite)]
        debris_objects = [obj for obj in self.space_objects if isinstance(obj, Debris)]

        if not os.path.exists("output"):
            os.makedirs("output")

        for i, sat in enumerate(satellites):
            print(f"Storing the satellite positions and miss distances data to {f'satellite_{i+1}.csv'}")
            with open(f"output/satellite_{i+1}.csv", "w", newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(["time", "x", "y", "z"] + [f"debris_{j+1}" for j in range(len(debris_objects))])
                for time_step in time_steps:
                    sat_distances = []
                    for deb in debris_objects:
                        sat_pos = sat.propagate(time_step * u.s).r.to(u.km).value
                        deb_pos = deb.propagate(time_step * u.s).r.to(u.km).value
                        distance = np.linalg.norm(sat_pos - deb_pos)
                        sat_distances.append(distance)
                    csv_writer.writerow([time_step, sat_pos[0], sat_pos[1], sat_pos[2]] + sat_distances)

class ConstellationVisualizer:
    def __init__(self, constellation):
        self.constellation = constellation
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.previous_positions = {i: [] for i in range(len(self.constellation.space_objects))}

    def plot_earth(self):
        """Plot Earth as a sphere."""
        earth_radius = Earth.R.to(u.km).value
        a, v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:25j]
        x = earth_radius * np.cos(a) * np.sin(v)
        y = earth_radius * np.sin(a) * np.sin(v)
        z = earth_radius * np.cos(v)
        self.ax.plot_surface(x, y, z, color='blue', alpha=0.4)

    def update(self, frame):
        """Update the constellation visualization for the given frame."""
        time = Time("J2000") + frame * u.s
        self.ax.clear()
        self.plot_earth()
        for i, space_object in enumerate(self.constellation.space_objects):
            propagated_orbit = space_object.propagate(frame * u.s)
            r = propagated_orbit.r.to(u.km).value
            marker_style = '2' if isinstance(space_object, Satellite) else 'o'
            s = '80' if isinstance(space_object, Satellite) else '20'
            label = f"Satellite {i+1}" if isinstance(space_object, Satellite) else f"Debris {i+1-8}"
            self.ax.scatter(r[0], r[1], r[2], marker=marker_style, label=label, s=int(s))
            self.ax.legend()

            # Store and plot previous positions
            self.previous_positions[i].append(r)
            if len(self.previous_positions[i]) > 1:
                prev_r = np.array(self.previous_positions[i])
                self.ax.plot(prev_r[:, 0], prev_r[:, 1], prev_r[:, 2], marker='.', linestyle='dotted', alpha=0.3, markersize=3)

        self.ax.set_xlabel("X (km)")
        self.ax.set_ylabel("Y (km)")
        self.ax.set_zlabel("Z (km)")
        self.ax.set_title(f"Constellation Propagation (t = {frame:.0f} s)")
        self.ax.set_xlim(-8000, 8000)
        self.ax.set_ylim(-8000, 8000)
        self.ax.set_zlim(-8000, 8000)

    def visualize(self, propagation_duration, time_steps):
        """Visualize the constellation propagation."""
        ani = FuncAnimation(self.fig, self.update, frames=time_steps, interval=1, blit=False)
        # ani.save('image/SCSA.gif')
        plt.show()

    def save_positions_and_distances(self, time_steps):
        """Save satellite positions and miss distances to separate CSV files."""
        self.constellation.save_positions_and_distances_to_csv(time_steps)


if __name__ == "__main__":
    # Define constellation parameters
    N = 8 # Number of satellites
    altitude = 700 * u.km
    num_planes = 4 # Number of orbital planes
    N_debris = 10
    debris_altitude = 700 * u.km
    debris_inclination_range = (95, 100) # inclination range in degrees

    # Create a constellation with satellites and debris
    constellation = Constellation()
    constellation.create_sso_constellation(N, altitude, num_planes)
    constellation.create_space_debris(N_debris, debris_altitude, debris_inclination_range)

    # Plot constellation
    propagation_duration = 1 * u.hour
    time_steps = np.linspace(0, propagation_duration.to(u.s).value, 360)

    visualizer = ConstellationVisualizer(constellation)
    visualizer.visualize(propagation_duration, time_steps)
    visualizer.save_positions_and_distances(time_steps)

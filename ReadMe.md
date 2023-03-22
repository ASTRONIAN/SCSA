# Space Constellation Simulation and Analysis

 This Python script simulates a constellation of satellites and space debris in Earth orbit. It propagates their orbits over time and calculates the miss distances between satellites and debris. The script also saves the positions of the satellites and miss distances as separate CSV files.
 
 ![SSO Constellation and Space Debris](https://github.com/ASTRONIAN/SCSA/blob/master/image/SCSA.gif)
 
## Dependencies
1. numpy
2. csv
3. os
4. astropy
5. poliastro
6. matplotlib
7. astropy

## SpaceObject Class
The SpaceObject class is a base class for all space objects, including satellites and debris. It initializes the object with its orbital parameters and provides a method to propagate the object's orbit over time.

## Satellite and Debris Classes
The Satellite and Debris classes inherit from the SpaceObject class. They define specific satellite and debris objects, respectively.

## Constellation Class
The Constellation class manages a collection of space objects. It provides methods to create a Sun-synchronous orbit (SSO) constellation, generate space debris, calculate miss distances between satellites and debris, and save satellite positions and miss distances to separate CSV files.

## ConstellationVisualizer Class
The ConstellationVisualizer class is responsible for visualizing the propagation of the space objects in the constellation. It provides methods to plot the Earth, update the positions of the space objects, visualize the constellation, and save satellite positions and miss distances to separate CSV files.

## Main Script
The main script defines the constellation parameters, creates a constellation with satellites and debris, and propagates the orbits over time. The script visualizes the constellation or saves the satellite positions and miss distances as separate CSV files.

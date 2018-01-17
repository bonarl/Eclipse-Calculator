# Eclipse-Calculator
Simulates the orbit of a satellite around mars and calculates the time the satellite spends in mars' shadow cone.

Satellite orbit is defined by orbital period, radius of periares and inclination. Can also define initial position
of mars by adjusting mars_0_angle, mars orbit is apprximated as circular. sim_time is running time of simulation
in hours, data_step the time between measurements.

Program outputs a graph of angle between mars-satellite vector and mars-sun vector (angle of satellite above
ecliptic plane) in blue, and red region at bottom of graph is the angle which the satellite would be in the
shadow cone, for the satellites altitude at time of measurement. The maximum eclipse duration, time and number
of eclipses are also output.

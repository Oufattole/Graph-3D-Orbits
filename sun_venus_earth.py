###BOOKLISTSTART1###
from amuse.lab import Particles, units 
import pandas as pd
import numpy as np

def new_system():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (855251, -804836, -3186) |units.km
    sun.velocity = (7.893, 11.894, 0.20642) | (units.m/units.s)

    venus = particles[1]
    venus.mass = 0.0025642 | units.MJupiter
    venus.radius = 3026.0 | units.km
    venus.position = (-0.3767, 0.60159, 0.03930) | units.AU
    venus.velocity = (-29.7725, -18.849, 0.795) | units.kms

    earth = particles[2]
    earth.mass = 1.0 | units.MEarth
    earth.radius = 1.0 | units.REarth
    earth.position = (-0.98561, 0.0762, -7.847e-5) | units.AU  
    earth.velocity = (-2.927, -29.803, -0.0005327) | units.kms

    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.lab import Hermite, nbody_system, Brutus
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(),
                                             particles[1].position.length())

    gravity = Hermite(convert_nbody)
    gravity.particles.add_particles(particles)
    sun = gravity.particles[0]
    venus = gravity.particles[1]
    earth = gravity.particles[2]
    def to_float(c):
        return [c[0].value_in(units.AU), c[1].value_in(units.AU), c[2].value_in(units.AU)]
    sun_pos = [to_float([sun.x,sun.y,sun.z])]
    earth_pos = [to_float([earth.x,earth.y,earth.z])]
    venus_pos = [to_float([venus.x,venus.y,venus.z])]
    start_array = np.array([sun_pos, earth_pos, venus_pos]).reshape([-1,3])
    start = pd.DataFrame(data=start_array, columns=["x", "y", "z"])
    
    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (1 | units.day))
        sun_pos.append(to_float([sun.x,sun.y,sun.z]))
        earth_pos.append(to_float([earth.x,earth.y,earth.z]))
        venus_pos.append(to_float([venus.x,venus.y,venus.z]))
    gravity.stop()
    return sun_pos, earth_pos, venus_pos, start
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(planet_1, planet_2, planet_3, start, output_filename):
    import matplotlib.pyplot as plot
    from mpl_toolkits.mplot3d import Axes3D

    x_label = 'x [au]'
    y_label = 'y [au]'
    z_label = 'z [au]'

    fig = plot.figure(figsize=(16, 12), dpi=80)

    colors = ['r','b','g']
    trajectories=[planet_1, planet_2, planet_3]
    ax = fig.add_subplot(1,1,1, projection='3d')
    max_range = 0
    for body, color in zip(trajectories, colors):
        body = np.array(body).reshape([-1,3])
        body = pd.DataFrame(data=body, columns=["x", "y", "z"])
        max_dim = max(max(body["x"]),max(body["y"]),max(body["z"]))
        if max_dim > max_range:
            max_range = max_dim
        ax.plot(body["x"], body["y"], body["z"], color)

    ax.scatter(start["x"], start["y"], start["z"], marker='o', s=200)

    ax.set_xlim([-max_range,max_range])    
    ax.set_ylim([-max_range,max_range])
    ax.set_zlim([-max_range,max_range])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    ax.legend()       

    save_file = 'SunVenusEarth.png'
    plot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    plot.show()
###BOOKLISTSTOP3###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-o", 
                      dest="output_filename", default ="SunVenusEarth",
                      help="output filename [%default]")
    return result
    
if __name__ in ('__main__','__plot__'):
    o, arguments  = new_option_parser().parse_args()

    particles = new_system()
    sun, earth, venus, start = integrate_solar_system(particles, 4 | units.yr)
    plot_track(sun, earth, venus, start, o.output_filename)
    
from amuse.lab import Particles, units 
import pandas as pd
import numpy as np

###BOOKLISTSTART1###
def new_system():
    stars = Particles(3)
    sun = stars[0]
    sun.mass = units.MSun(1.0)
    sun.position = units.m(np.array((0.0,0.0,0.0)))
    sun.velocity = units.ms(np.array((0.0,0.0,0.0)))
    sun.radius = units.RSun(1.0)

    earth = stars[1]
    earth.mass = units.kg(5.9736e24)
    earth.radius = units.km(6371) 
    earth.position = units.km(np.array((149.5e6,0.0,0.0)))
    earth.velocity = units.ms(np.array((0.0,29800,0.0)))

    moon = stars[2]
    moon.mass = units.kg(7.3477e22 )
    moon.radius = units.km(1737.10) 
    moon.position = units.km(np.array((149.5e6 + 384399.0 ,0.0,0.0)))
    moon.velocity = ([0.0,1.022,0] | units.km/units.s) + earth.velocity    
    return stars
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.lab import Hermite, nbody_system
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
    
    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (1 | units.day))
        sun_pos.append(to_float([sun.x,sun.y,sun.z]))
        earth_pos.append(to_float([earth.x,earth.y,earth.z]))
        venus_pos.append(to_float([venus.x,venus.y,venus.z]))
    gravity.stop()
    return sun_pos, earth_pos, venus_pos
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
    start_array = np.array([earth_pos, venus_pos]).reshape([-1,3])
    start = pd.DataFrame(data=start_array, columns=["x", "y", "z"])
    
    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (.01 | units.day))
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
    planet_2 = pd.DataFrame(data=np.array(planet_2).reshape([-1,3]), columns=["x", "y", "z"])
    planet_3 = pd.DataFrame(data=np.array(planet_3).reshape([-1,3]), columns=["x", "y", "z"])
    trajectories=[planet_2, planet_3]
    ax = fig.add_subplot(1,1,1, projection='3d')
    max_rangex = max(max(planet_2["x"]),max(planet_3['x']))
    min_rangex = min(min(planet_2["x"]),min(planet_3['x']))
    max_rangey = max(max(planet_2["y"]),max(planet_3['y'])) 
    min_rangey = min(min(planet_2["y"]),min(planet_3['y']))
    max_rangez = max(max(planet_2["z"]),max(planet_3['z'])) 
    min_rangez = min(min(planet_2["z"]),min(planet_3['z']))
    for body, color in zip(trajectories, colors):

        ax.plot(body["x"], body["y"], body["z"], color)

    ax.scatter(start["x"], start["y"], start["z"], marker='o', s=200)

    ax.set_xlim([min_rangex,max_rangex]) 
    ax.set_ylim([min_rangey,max_rangey])
    ax.set_zlim([min_rangez,max_rangez])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    ax.legend()       

    save_file = 'EarthMoonSun.png'
    plot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    plot.show()
###BOOKLISTSTOP3###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-o", 
                      dest="output_filename", default ="EarthMoonSun",
                      help="output filename [%default]")
    return result
    
if __name__ in ('__main__','__plot__'):
    o, arguments  = new_option_parser().parse_args()

    particles = new_system()
    sun, earth, venus, start = integrate_solar_system(particles, .1 | units.yr)
    plot_track(sun, earth, venus, start, o.output_filename)
    
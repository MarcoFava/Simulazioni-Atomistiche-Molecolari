import numpy

class Particle:

    # def __init__(self, position, velocity=None, species='Ar', mass=1.0):
    #     self.position = numpy.asarray(position)
    #     if velocity is not None:
    #         self.velocity = numpy.asarray(velocity)
    #     else:
    #         self.velocity = numpy.zeros(len(position))
    #     self.species = species
    #     self.mass = mass

        

    def __init__(self, position, velocity=[0.0, 0., 0.], species='Ar', mass=1.0):
        self.position = numpy.asarray(position)
        self.velocity = velocity
        # if velocity is not None:
        #     self.velocity = numpy.asarray(velocity)
        # else:
        #     self.velocity = numpy.zeros(len(position))
        self.species = species
        self.mass = mass
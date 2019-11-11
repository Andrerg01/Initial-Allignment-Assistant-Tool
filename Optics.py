import numpy as np

#Returns the norm/magnitude of a N-Dimensional input vector in cartesian coordinates
def norm(vector):
    v = np.array(vector)
    return np.sqrt(v.dot(v))

#Returns a normalized version of the N-Dimensional input vector in cartesian coordinates
def normalize(vector):
    v = np.array(vector)
    return v/norm(v)

#Returns an array with the values of r, th, and ph for the 3-Dimensional input vector in cartesian coordinates (r = norm of vector, th = polar angle (going from -pi2/ to pi/2, ph = azymuthal angle)
def toSpherical(vector):
    r = norm(vector)
    th = np.arctan2(vector[2], np.sqrt(vector[0]**2+vector[1]**2))
    ph = np.arctan2(vector[1], vector[0])      
    return np.array([r,th,ph])

#Return a 3-Dimensional vector in cartesian coordinates given an array with the values of r, th and ph (r = norm of vector, th = polar angle (from -pi/2 to pi/2), ph = azymuthal angle)
def toCartesian(sphericalArray):
    r = sphericalArray[0]
    th = sphericalArray[1]
    ph = sphericalArray[2]
    x = r*np.cos(th)*np.cos(ph)
    y = r*np.cos(th)*np.sin(ph)
    z = r*np.sin(th)
    return np.array([x,y,z])

#Returns the angle between two vectors in cartesian coordinates.
def angleBetweenVectors(vector1, vector2):
    return np.arccos(vector1.dot(vector2)/(norm(vector1)*norm(vector2)))

#Returns the critical angle for refraction.
def criticalAngle(indexOfRefraction):
    return arcsin(1.0/indexOfRefraction)

#Changes a vector's "Yaw" by changing its azymuthal angle, returns a rotated version of the vector
def rotateYaw(vector, yaw):
    vectorSph = toSpherical(vector)
    vectorSph[2] = vectorSph[2] + yaw
    return toCartesian(vectorSph)

#Changes a vector's "Pitch" by changing its polar angle, returns a rotated version of the vector
def rotatePitch(vector, pitch):
    vectorSph = toSpherical(vector)
    vectorSph[1] = vectorSph[1] + pitch
    return toCartesian(vectorSph)

#Changes a vector pitch and yaw by changing its polar and azymuthal angles, the order does no matter
def rotatePitchYaw(vector, pitch, yaw):
    return rotatePitch(rotateYaw(vector, yaw), pitch)

#Returns True of False if the beam will eventually collide with the optical element
#This works for both lenses and mirrors    
def intersectionBetweenLineAndSphere(directionOfLine, pointInLine, centerOfSphere, radiusOfSphere):
    #Applying formula for intersection between line and sphere.
    l = directionOfLine
    o = pointInLine
    c = centerOfSphere
    r = radiusOfSphere
    #This is the part of the distance formula that can give a complex value, meaning there is no collision.
    if (l.dot(o-c))**2 - (norm(o-c)**2-r**2) > 0:
        #The two points where the line collides the sphere are.
        d1 = -(l.dot(o-c)) + np.sqrt((l.dot(o-c))**2 - (norm(o-c)**2-r**2))
        d2 = -(l.dot(o-c)) - np.sqrt((l.dot(o-c))**2 - (norm(o-c)**2-r**2))
        if r >= 0:
            #pick the greatest of the distances, which will also be the one that intersects the sphere form the inside, for a concave mirror (r>=0).
            d = max([d1,d2])
        else:
            #pick the smallest of the distances, which will also be the one that intersects the sphere form the outside, for a convex mirror (r<0).
            d = min([d1,d2])
        #The intersection point will be the line of the beam after travelling a distance d
        intersection = o + d*l
        return intersection
    else:
        return None

def intersectionBetweenLineAndPlane(directionOfLine, pointInLine, pointInPlane, normalOfPlane):
    p = pointInPlane
    n = normalOfPlane
    l = directionOfLine
    o = pointInLine
    if l.dot(n) != 0:
        d = (p - o).dot(n)/(l.dot(n))
        intersection = o + d*l
    elif (p-o).dot(n) == 0:
        intersection = o
    else:
        intersection = None
    return intersection
        
#Class used for the propagation of a Gaussian Beam    
class Beam:
    #Initialization of the class    
    def __init__(self, radiusOfCurvature, width, direction, position = [0,0,0], wavelength = 1064.0E-9, indexOfRefraction = 1, verbose = True):
        #Current radius of curvature of the beam
        self.radiusOfCurvature = radiusOfCurvature
        #Initial width (sigma) of the beam
        self.width = width
        #Wavelength of the beam
        self.wavelength = wavelength
        #Index of refraction of the medium
        self.indexOfRefraction = indexOfRefraction
        #Current position of the beam (given as an array or a np.array of 3 values, [x,y,z]), code converts to a np.array
        self.position = np.array(position)
        #Current direction of the beam (given as an arrat or a np.array of 3 values[x,y,z]), code converts to a np.arrat and normalizes it. (essentially k-vector)
        self.direction = normalize(direction)
        #Toggles diologs when interactions happens
        self.verbose = verbose
    #Returns the parameter q such that 1/q = 1/r - i*lambda/(pi*n*w^2)
    def q(self):
        #numpy used j as the imaginary unit
        return 1.0/(1.0/self.radiusOfCurvature - 1.0j*self.wavelength/(np.pi*self.indexOfRefraction*self.width**2.0))

    #Returns the radius of curvature form a given q    
    def radiusFrom_q(self,q):
        return ((q**(-1.0)).real)**(-1.0)
    #Returns the width form a given q
    def widthFrom_q(self,q):
        return (((q**(-1.0)).imag)*(-np.pi*self.indexOfRefraction/self.wavelength))**(-0.5)
  
    #Returns a copy of the beam state as it is when the function is called    
    def copy(self):
        return Beam(radiusOfCurvature = self.radiusOfCurvature, width = self.width, direction = self.direction, position = self.position, wavelength = self.wavelength, indexOfRefraction = self.indexOfRefraction, verbose = self.verbose)
    
    #Propagates the beam for a given distance. Will update the value of the radius of cruvature, width, and position. It will propagate the beam in the direction that it has. Does not return anything
    def propagate(self, distance):
        if self.verbose: 
            print("Propagating a distance of " + str(distance))
        q = self.q()
        q = q + distance
        self.radiusOfCurvature = self.radiusFrom_q(q)
        self.width = self.widthFrom_q(q)
        self.position = self.position + distance*self.direction
        
    #Changes the direction of the beam after reflecting on a plane perpendicular to the normal being inputed
    def reflect(self, normal):
        self.direction = self.direction - 2*(self.direction.dot(normal))*normal
    
    def refract(self, element):
        n1 = 1.0
        n2 = 1.0*element.indexOfRefraction
        k = n1/n2
        normal = normalize(element.center1() - self.position)
        thI = angleBetweenVectors(-1*self.direction, normal)
        self.direction = normalize(k*self.direction + (k*np.cos(thI)-np.sqrt(1.0-(k**2)*(1.0-np.cos(thI)**2.0)))*normal)
        self.indexOfRefraction = n2
        q = self.q()
        q = 1.0/(1.0/q - (n2-n1)/element.radiusOfCurvature)
        self.radiusOfCurvature = self.radiusFrom_q(q)
        self.width = self.widthFrom_q(q)
        distanceToOtherSide = norm(self.position - intersectionBetweenLineAndSphere(self.direction, self.position, element.center2(), -element.radiusOfCurvature))
        self.propagate(distanceToOtherSide)
        n1 = 1.0*element.indexOfRefraction
        n2 = 1.0
        k = n1/n2
        normal = -normalize(element.center2() - self.position)
        self.direction = normalize(k*self.direction + (k*np.cos(thI)-np.sqrt(1.0-(k**2)*(1.0-np.cos(thI)**2.0)))*normal)
        self.indexOfRefraction = n2
        q = self.q()
        q = 1.0/(1.0/q - (n2-n1)/element.radiusOfCurvature)
        self.radiusOfCurvature = self.radiusFrom_q(q)
        self.width = self.widthFrom_q(q)
        
    #Returns True of False if the beam will eventually collide with the optical element
    #This works for both lenses and mirrors    
    def collisionQ(self, element):
        if isinstance(element, Mirror) or isinstance(element, Lens):
            #Applying formula for intersection between line and sphere.
            intersection = intersectionBetweenLineAndSphere(self.direction, self.position, element.center1(), element.radiusOfCurvature)
            #Intersection function will return None if there is no internsection
        elif isinstance(element, FlatMirror) or isinstance(element, InfinitePlane):
            intersection = intersectionBetweenLineAndPlane(self.direction, self.position, element.vertex(), element.normal())
        
        if intersection is None:
            return False
        #Makes sure the intersection point is withing the diameter of the element (not considering spherical caps here), and makes sure it is in front of the beam's path, not behind.
        if isinstance(element, InfinitePlane):
            return True
        if norm(intersection - element.vertex1()) <= element.diameter/2.0 and self.direction.dot(intersection - self.position) >= 0:
            return True
            
        return False
        
    #Returns True of False if the beam will eventually collide with the optical element
    #This works for both lenses and mirrors    
    def collisionPoint(self, element):
        if isinstance(element, Mirror) or isinstance(element, Lens):
            #Applying formula for intersection between line and sphere.
            intersection = intersectionBetweenLineAndSphere(self.direction, self.position, element.center1(), element.radiusOfCurvature)
            #Intersection function will return None if there is no internsection
        elif isinstance(element, FlatMirror) or isinstance(element, InfinitePlane) or isinstance(element, Aperture):
            intersection = intersectionBetweenLineAndPlane(self.direction, self.position, element.vertex(), element.normal())
        if intersection is None:
            return False
        #Makes sure the intersection point is withing the diameter of the element (not considering spherical caps here), and makes sure it is in front of the beam's path, not behind.
        if isinstance(element, InfinitePlane) or norm(intersection - element.vertex1()) <= element.diameter/2 and self.direction.dot(intersection - self.position) >= 0:
            return intersection
        else:
            return None

    def clippingQ(self, element):
        if norm(element.vertex1() - self.position) + self.width > element.diameter/2.0:
            return True
        else:
            return False
    
    def interact(self, element):
        if self.collisionQ(element):
            distanceToElement = norm(self.position - self.collisionPoint(element))
            self.propagate(distanceToElement)
            if isinstance(element, Mirror):
                normalToCollisionPoint = normalize(element.center1() - self.position)
                self.reflect(normalToCollisionPoint)
                q = self.q()
                q = 1.0/(1.0/q - 2.0/element.radiusOfCurvature)
                self.radiusOfCurvature = self.radiusFrom_q(q)
                self.width = self.widthFrom_q(q)
            elif isinstance(element, Lens):
                self.refract(element)
            elif isinstance(element, FlatMirror):
                self.reflect(element.normal())
            elif isinstance(element, InfinitePlane):
                self.reflect(element.normal())
        else:
            if self.verbose: print("Cannot interact with element " + str(element.ID) + ", no collision detected!")
            self.interact(InfinitePlane(element.ID, element.vertex1(), element.yaw, element.pitch))
                
    def track(self, elements, step):
        beam = self.copy()
        beamStates = np.array([beam.copy()])
        for element in elements:
            steps = np.int(norm(beam.position-beam.collisionPoint(element))/step)-1
            for i in range(steps):
                beam.propagate(step)
                beamStates = np.append(beamStates, beam.copy())
            beam.interact(element)
            beamStates = np.append(beamStates, beam.copy())
        for i in range(30):
            beam.propagate(step)
            beamStates = np.append(beamStates, beam.copy())
        return beamStates
    
    def calculateStates(self, elements):
        beam = self.copy()
        beamStates = {'Source':beam.copy()}
        for element in elements:
            beam.interact(element)
            beamStates[element.ID] = beam.copy()
        return beamStates
    
    def calculateFlags(self, elements):
        beam = self.copy()
        flags = []
        for element in elements:
            if beam.collisionQ(element):
                beam.interact(element)
                if beam.clippingQ(element):
                    flags.append("Warning: Beam clipping with element " + str(element.ID))
            else:
                flags.append("Error!: Beam not intersecting with element " + str(element.ID))
        return flags
    
    #Nicely prints all the attributes of the beam at the moment the function is called
    def __str__(self):
        return \
        "Radius Of Curvature : " + str(self.radiusOfCurvature) + "\n" + \
        "Width : " + str(self.width) + "\n" + \
        "Parameter q : " + str(self.q()) + "\n" + \
        "Wavelength : " + str(self.wavelength) + "\n" + \
        "Index of Refraction : " + str(self.indexOfRefraction) + "\n" + \
        "Position : " + str(self.position) + "\n" + \
        "Direction : " + str(self.direction) + "\n"    
    
#Class used for the interaction of the gaussian beam with mirror elements"""    
class Mirror:
    #Initialization for the class
    def __init__(self, ID, radiusOfCurvature, positionOfCM, parameter_d, yaw, pitch, diameter):
        #ID of the mirror (for identification).
        self.ID = ID
        #Radius of curvatire of the mirror (2F).
        self.radiusOfCurvature = radiusOfCurvature
        #Position of the center of mass for which the mirror will rotate about for pitch and yaw adjustments.
        self.positionOfCM = np.array(positionOfCM)
        #Distance from the center of mass to the vertex of the mirror.
        self.parameter_d = parameter_d
        #Pitch of the mirror from (1,0,0) (pitch = polar rotation).
        self.pitch = pitch
        #Yaw of the mirror from (1,0,0) (yaw = axymuthal rotation)
        self.yaw = yaw
        #Diameter of the mirror
        self.diameter = diameter
        
    #Return the position of the center of the sphere of which the mirror is a cap of (think of the mirror as a section of big sphere, this is the center of that sphere), the '1' is so the same function can also be called for lenses without having to test for the type beforehand.
    def center1(self):
        return self.vertex1() + self.radiusOfCurvature*self.normal()
    
    def center(self):
        return self.center1()
    #Return the normal of the mirror given the pitch and yaw from (1,0,0)
    def normal(self):
        return rotatePitchYaw([1,0,0], self.pitch, self.yaw)

    #Returns the position of the actual vertex of the mirror, taking into account the deviation of it from the center of mass and yaw+pitch rotations, the '1' is so the same function can also be called for lenses without having to test for the type beforehand
    def vertex1(self):
        return self.positionOfCM + self.parameter_d*self.normal()
    
    def vertex(self):
        return self.vertex1()
    
    #Returns a copy of the mirror state as it is when the function is called
    def copy(self):
        return Mirror(ID = self.ID, radiusOfCurvature = self.radiusOfCurvature, positionOfCM = self.positionOfCM, parameter_d = self.parameter_d, yaw = self.yaw, pitch = self.pitch, diameter = self.diameter)
    
    #Nicely prints all the attributes of the mirror at the moment the function is called
    def __str__(self):
        return \
        "ID : " + str(self.ID) + "\n" + \
        "Radius Of Curvature : " + str(self.radiusOfCurvature) + "\n" + \
        "Position of Center of Mass : " + str(self.positionOfCM) + "\n" + \
        "Parameter_d : " + str(self.parameter_d) + "\n" + \
        "yaw : " + str(self.yaw) + "\n" + \
        "Pitch : " + str(self.pitch) + "\n" + \
        "Normal : " + str(self.normal()) + "\n" + \
        "Diameter : " + str(self.diameter) + "\n"    

class FlatMirror:
    def __init__(self, ID, positionOfCM, parameter_d, yaw, pitch, diameter):
        #ID of the flat mirror
        self.ID = ID
        #Position of the center of mass for which the flat mirror will rotate about for pitch and yaw adjustments
        self.positionOfCM = positionOfCM
        #Distance between the vertex (surface) of the flat mirror and the center of mass
        self.parameter_d = parameter_d
        #Pitch of the flat mirror from (1,0,0) (pitch = polar rotation).
        self.pitch = pitch
        #Yaw of the flat mirror from (1,0,0) (yaw = axymuthal rotation)
        self.yaw = yaw
        #Diameter of the Lens
        self.diameter = diameter
    def normal(self):
        return rotatePitchYaw([1,0,0], self.pitch, self.yaw)
    def vertex1(self):
        return self.positionOfCM + self.parameter_d*self.normal()
    def vertex(self):
        return self.vertex1()
    def copy(self):
        return FlatMirror(ID = self.ID, positionOfCM = self.positionOfCM, parameter_d = self.parameter_d, yaw = self.yaw, pitch = self.pitch, diameter = self.diameter)
    def __str__(self):
        return \
        "ID : " + str(self.ID) + "\n" + \
        "Position of Center of Mass : " + str(self.positionOfCM) + "\n" + \
        "Parameter_d : " + str(self.parameter_d) + "\n" + \
        "yaw : " + str(self.yaw) + "\n" + \
        "Pitch : " + str(self.pitch) + "\n" + \
        "Normal : " + str(self.normal()) + "\n" + \
        "Diameter : " + str(self.diameter) + "\n"
    
class Lens:
    def __init__(self, ID, radiusOfCurvature, positionOfCM, parameter_d, yaw, pitch, diameter, indexOfRefraction):
        #ID of the lens (for identification)
        self.ID = ID
        #Radius of curvature of the lens, assuming spherical, so 2F for all intents and purposes, R > 0 for concave, R < 0 for convex.
        self.radiusOfCurvature = radiusOfCurvature
        #Position of the center of mass for which the lens will rotate about for pitch and yaw adjustments.
        self.positionOfCM = np.array(positionOfCM)
        #Distance between the vertex of the lens and the center of mass
        self.parameter_d = parameter_d
        #Pitch of the lens from (1,0,0) (pitch = polar rotation).
        self.pitch = pitch
        #Yaw of the lens from (1,0,0) (yaw = axymuthal rotation)
        self.yaw = yaw
        #Diameter of the Lens
        self.diameter = diameter
        #Index of refraction inside the lens.
        self.indexOfRefraction = indexOfRefraction
    
    #Returns the position of the center of the sphere of which the side of the lens where the normal is calculated is a cap of.
    def center1(self):
        return self.vertex1() + self.radiusOfCurvature*self.normal()
    #Returns the position of the center of the sphere of which the opposide side of the lens where the vertex is NOT located is a cap of.
    def center2(self):
        return self.vertex2() - self.radiusOfCurvature*self.normal()
    
    #Return the normal of the lens given the pitch and yaw from (1,0,0)
    def normal(self):
        return rotatePitchYaw([1,0,0], self.pitch, self.yaw)
    
    #Returns the position of the vertex of the lens where is normal is calculated from.
    def vertex1(self):
        return self.positionOfCM + self.parameter_d*self.normal()
    
    #Returns the position of the vertex of the lens where is normal NOR is calculated from. (opposite to vertex1)
    def vertex2(self):
        return self.positionOfCM - self.parameter_d*self.normal()
    
    #Returns a copy of the mirror state as it is when the function is called
    def copy(self):
        return Lens(ID = self.ID, radiusOfCurvature = self.radiusOfCurvature, positionOfCM = self.positionOfCM, parameter_d = self.parameter_d, yaw = self.yaw, pitch = self.pitch, diameter = self.diameter, indexOfRefraction = self.indexOfRefraction)
 
    #Nicely prints all the attributes of the mirror at the moment the function is called
    def __str__(self):
        return \
        "ID : " + str(self.ID) + "\n" + \
        "Radius Of Curvature : " + str(self.radiusOfCurvature) + "\n" + \
        "Position of Center of Mass : " + str(self.positionOfCM) + "\n" + \
        "Parameter_d : " + str(self.parameter_d) + "\n" + \
        "yaw : " + str(self.yaw) + "\n" + \
        "Pitch : " + str(self.pitch) + "\n" + \
        "Normal : " + str(self.normal()) + "\n" + \
        "Diameter : " + str(self.diameter) + "\n" + \
        "Index of Refraction : " + str(self.indexOfRefraction) + "\n"

class Aperture:
    def __init__(self, ID, positionOfCM, yaw, pitch, diameter):
        #ID of the flat mirror
        self.ID = ID
        #Position of the center of mass for which the flat mirror will rotate about for pitch and yaw adjustments
        self.positionOfCM = positionOfCM
        #Pitch of the flat mirror from (1,0,0) (pitch = polar rotation).
        self.pitch = pitch
        #Yaw of the flat mirror from (1,0,0) (yaw = axymuthal rotation)
        self.yaw = yaw
        #Diameter of the Lens
        self.diameter = diameter
    def normal(self):
        return rotatePitchYaw([1,0,0], self.pitch, self.yaw)
    def vertex1(self):
        return self.positionOfCM
    def vertex(self):
        return self.vertex1()
    def copy(self):
        return Aperture(ID = self.ID, positionOfCM = self.positionOfCM, yaw = self.yaw, pitch = self.pitch, diameter = self.diameter)
    def __str__(self):
        return \
        "ID : " + str(self.ID) + "\n" + \
        "Position of Center of Mass : " + str(self.positionOfCM) + "\n" + \
        "yaw : " + str(self.yaw) + "\n" + \
        "Pitch : " + str(self.pitch) + "\n" + \
        "Normal : " + str(self.normal()) + "\n" + \
        "Diameter : " + str(self.diameter) + "\n"
    
class InfinitePlane:
    def __init__(self, ID, positionOfCM, yaw, pitch):\
        #ID of the flat mirror
        self.ID = ID
        #Position of the center of mass for which the flat mirror will rotate about for pitch and yaw adjustments
        self.positionOfCM = positionOfCM
        #Pitch of the flat mirror from (1,0,0) (pitch = polar rotation).
        self.pitch = pitch
        #Yaw of the flat mirror from (1,0,0) (yaw = axymuthal rotation)
        self.yaw = yaw
    def normal(self):
        return rotatePitchYaw([1,0,0], self.pitch, self.yaw)
    def vertex1(self):
        return self.positionOfCM
    def vertex(self):
        return self.vertex1()
    def copy(self):
        return InfinitePlane(ID = self.ID, positionOfCM = self.positionOfCM, yaw = self.yaw, pitch = self.pitch)
    def __str__(self):
        return \
        "ID : " + str(self.ID) + "\n" + \
        "Position of Center of Mass : " + str(self.positionOfCM) + "\n" + \
        "yaw : " + str(self.yaw) + "\n" + \
        "Pitch : " + str(self.pitch) + "\n" + \
        "Normal : " + str(self.normal()) + "\n"
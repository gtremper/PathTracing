/* Class description for Ray */
#ifndef RAY_H
#define RAY_H

class Ray {
public:
	Ray(){};
	Ray(const vec3&,const vec3&);
	
	vec3 getPoint(float);
	vec3 origin;
	vec3 direction;	
};

struct pixel {
	float red;
	float green;
	float blue;
};

#endif
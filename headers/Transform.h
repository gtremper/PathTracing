// Transform header file to define the interface. 
// The class is all static for simplicity
// You need to implement left, up and lookAt
// Rotate is a helper function

// Include the helper glm library, including matrix transform extensions

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// glm provides vector, matrix classes like glsl
// Typedefs to make code more readable 

typedef glm::mat2 mat2 ;
typedef glm::mat3 mat3 ;
typedef glm::mat4 mat4 ; 
typedef glm::vec2 vec2 ;
typedef glm::vec3 vec3 ; 
typedef glm::vec4 vec4 ; 

const float pi = 3.14159265358979 ; // For portability across platforms

class Transform	 
{
public:
	Transform();
	virtual ~Transform();
	static mat4 rotate(float degrees, const vec3& axis) ;
	static mat4 scale(float sx, float sy, float sz) ; 
	static mat4 translate(float tx, float ty, float tz);
};

#endif


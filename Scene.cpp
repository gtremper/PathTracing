/* Deals with the rendered scene */
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <string>
#include <vector>
#include "float.h"

#include "FreeImage.h"
#include "Shapes.h"
#include "Intersection.h"
#include "Light.h"
#include "Scene.h"

using namespace std;


/***  SCENE  ***/

// sets default values
Scene::Scene(char* file) {
	filename = "OUTPUT.png";
	maxdepth = 5; 
	ambient = vec3(0.2,0.2,0.2);
	diffuse = vec3(0,0,0);
	specular = vec3(0,0,0);
	shininess = 0;
	emission = vec3(0,0,0);
	indexofrefraction = 1;
	refractivity = 0;
	antialias = 2;
	shadowrays = 1;
	lightradius = 0;
	
	constant = 1;
	linear = 0;
	quadratic = 0;
	vec3 maxvec = vec3(DBL_MAX,DBL_MAX,DBL_MAX);
	vec3 minvec = vec3(DBL_MIN,DBL_MIN,DBL_MIN);
	sceneAABB = AABB(maxvec, minvec);
	parse(file);
	clog << "Constructing KDTree... ";
	KDTree = new TreeNode(objects,sceneAABB);
	clog << "done"<<endl;
}

Scene::~Scene() {
	vector<Shape*>::iterator it;
	for(it=objects.begin(); it!=objects.end(); it++){
		delete *it;
	}
	vector<Light*>::iterator l;
	for(l=lights.begin(); l!=lights.end(); l++){
		delete *l;
	}
}
 
Ray Scene::castEyeRay(double i, double j){
	double alpha = (2.0*i-width)/width;
	alpha *=tan(fovx/2.0);
	double beta = (2.0*j-height)/height;
	beta *= tan(fovy/2.0);
	
	Ray ray(eye,glm::normalize(alpha*u + beta*v - w));
	return ray;
}

/***********************************************
** Parsing input and setting global variables **
***********************************************/

void Scene::setCoordinateFrame(vec3& lookat, vec3& up){
	w = glm::normalize(eye - lookat);
	u = glm::normalize(glm::cross(up, w));
	v = glm::cross(w,u);
}

void Scene::updateAABB(vec3& point){
	sceneAABB.aabbmin[0] = min(sceneAABB.aabbmin[0],point[0]);
	sceneAABB.aabbmax[0] = max(sceneAABB.aabbmax[0],point[0]);
	sceneAABB.aabbmin[1] = min(sceneAABB.aabbmin[1],point[1]);
	sceneAABB.aabbmax[1] = max(sceneAABB.aabbmax[1],point[1]);
	sceneAABB.aabbmin[2] = min(sceneAABB.aabbmin[2],point[2]);
	sceneAABB.aabbmax[2] = max(sceneAABB.aabbmax[2],point[2]);
}

void Scene::parseLine(string l, stack<mat4>& mv, vector<vec3>& verts, 
				vector<vec3>& normverts, vector<vec3>& norms){
					
	stringstream line(l);
	string cmd;
	line >> cmd;
	if (cmd[0]=='#' || cmd=="") { //comment or blank line
		return;
	} else if (cmd == "size") {
		line >> width >> height;
	} else if (cmd == "maxdepth") {
		line >> maxdepth;
	} else if (cmd == "output") {
		line >> filename;
	} else if (cmd == "camera") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		eye = vec3(arg1,arg2,arg3);
		
		line >> arg1 >> arg2 >> arg3;
		vec3 lookat = vec3(arg1,arg2,arg3);
		
		line >> arg1 >> arg2 >> arg3;
		vec3 up = vec3(arg1,arg2,arg3);
		setCoordinateFrame(lookat,up);
		line >> fovy;
		fovy*=pi/180.0;
		double d = height/(2.0 * tan(fovy*0.5));
		fovx = 2. * atan(width/(2.0*d));
	} else if (cmd == "sphere") {
		double arg1, arg2, arg3, arg4;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		line >> arg4;
		mat4 trans = mv.top();
		trans *= Transform::translate(arg1,arg2,arg3);
		trans *= Transform::scale(arg4,arg4,arg4);
		Sphere* s = new Sphere(trans);
		sceneAABB.aabbmin[0] = min(s->aabb.aabbmin[0],sceneAABB.aabbmin[0]);
		sceneAABB.aabbmax[0] = max(s->aabb.aabbmax[0],sceneAABB.aabbmax[0]);
		sceneAABB.aabbmin[1] = min(s->aabb.aabbmin[1],sceneAABB.aabbmin[1]);
		sceneAABB.aabbmax[1] = max(s->aabb.aabbmax[1],sceneAABB.aabbmax[1]);
		sceneAABB.aabbmin[2] = min(s->aabb.aabbmin[2],sceneAABB.aabbmin[2]);
		sceneAABB.aabbmax[2] = max(s->aabb.aabbmax[2],sceneAABB.aabbmax[2]);
		s->ambient = ambient;
		s->diffuse = diffuse;
		s->specular = specular;
		s->shininess = shininess;
		s->emission = emission;	
		s->emission = emission;
		s->indexofrefraction = indexofrefraction;
		s->refractivity = refractivity;
		objects.push_back(s);
	} else if (cmd == "maxverts") {
		int maxverts;
		line >> maxverts;
		verts.reserve(maxverts);
	} else if (cmd == "maxvertnorms") {
		int maxvertnorms;
		line >> maxvertnorms;
		verts.reserve(maxvertnorms);
		verts.reserve(maxvertnorms);
	} else if (cmd == "vertex") {
		double arg1, arg2, arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		vec3 v(arg1,arg2,arg3);
		verts.push_back(v);
	} else if (cmd == "vertexnormal") {
		double arg1, arg2, arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		vec3 v(arg1,arg2,arg3);
		normverts.push_back(v);
		line >> arg1;
		line >> arg2;
		line >> arg3;
		vec3 n(arg1,arg2,arg3);
		norms.push_back(n);
	} else if (cmd == "tri") {
		int a1, a2, a3;
		line >> a1 >> a2 >> a3;
		mat4 top =  mv.top();
		vec3 v1 = vec3(top * vec4(verts[a1],1));
		vec3 v2 = vec3(top * vec4(verts[a2],1));
		vec3 v3 = vec3(top * vec4(verts[a3],1));
		Triangle* t = new Triangle(v1,v2,v3);
		updateAABB(v1);
		updateAABB(v2);
		updateAABB(v3);
		t->ambient = ambient;
		t->diffuse = diffuse;
		t->specular = specular;
		t->shininess = shininess;
		t->emission = emission;
		t->indexofrefraction = indexofrefraction;
		t->refractivity = refractivity;
		objects.push_back(t);
	} else if(cmd == "trinormal") {
		int a1,a2,a3;
		line >> a1 >> a2 >> a3;
		mat4 top =  mv.top();
		vec3 v1 = vec3(top * vec4(normverts[a1],1));
		vec3 v2 = vec3(top * vec4(normverts[a2],1));
		vec3 v3 = vec3(top * vec4(normverts[a3],1));
		top = glm::transpose(glm::inverse(top));
		vec3 n1 = vec3(top * vec4(norms[a1],0));
		vec3 n2 = vec3(top * vec4(norms[a2],0));
		vec3 n3 = vec3(top * vec4(norms[a3],0));
		NormTriangle* t = new NormTriangle(v1,v2,v3,n1,n2,n3);
		t->ambient = ambient;
		t->diffuse = diffuse;
		t->specular = specular;
		t->shininess = shininess;
		t->emission = emission;
		t->indexofrefraction = indexofrefraction;
		t->refractivity = refractivity;
		objects.push_back(t);
	} else if(cmd == "translate") {
		double arg1,arg2,arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		mv.top() *= Transform::translate(arg1, arg2, arg3);
	} else if(cmd == "rotate") {
		double arg1,arg2,arg3,arg4;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		line >> arg4;
		mv.top() *= Transform::rotate(arg4,vec3(arg1,arg2,arg3));
	} else if (cmd=="scale") {
		double arg1,arg2,arg3;
		line >> arg1;
		line >> arg2;
		line >> arg3;
		mv.top() *= Transform::scale(arg1, arg2, arg3);
	} else if (cmd == "pushTransform") {
		mv.push(mv.top());
	} else if (cmd == "popTransform"){
		mv.pop();
	} else if (cmd == "directional") {
		double x,y,z,r,g,b;
		line >> x >> y >> z >> r >> g >> b;
		vec3 dir = vec3(mv.top()*vec4(x,y,z,0.0));
		DirectionalLight* light = new DirectionalLight(vec3(r,g,b),dir);
		lights.push_back(light);
	} else if (cmd == "point") {
		double x,y,z,r,g,b;
		line >> x >> y >> z >> r >> g >> b;
		vec3 point = vec3(mv.top()*vec4(x,y,z,1.0));
		PointLight* light = new PointLight(vec3(r,g,b),vec3(x,y,z),constant,linear,quadratic);
		light->shadowrays = shadowrays;
		light->lightradius = lightradius;
		lights.push_back(light);
	} else if (cmd == "attenuation") {
		line >> constant >> linear >> quadratic;
	} else if (cmd == "ambient") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		ambient = vec3(arg1,arg2,arg3);
	} else if (cmd == "diffuse") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		diffuse = vec3(arg1,arg2,arg3);
	} else if (cmd == "specular") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		specular = vec3(arg1,arg2,arg3);
	} else if (cmd == "shininess") {
		line >> shininess;
	} else if (cmd == "emission") {
		double arg1, arg2, arg3;
		line >> arg1 >> arg2 >> arg3;
		emission = vec3(arg1,arg2,arg3);
	} else if (cmd == "indexofrefraction") {
		line >> indexofrefraction;
	} else if (cmd == "refractivity") {
		line >> refractivity;
	} else if (cmd == "antialias") {
		line >> antialias;
	} else if (cmd == "shadowrays") {
		line >> shadowrays;
	} else if (cmd == "lightradius") {
		line >> lightradius;
	}
	//cout << cmd << endl;
}

void Scene::parse(char* filename) {
	ifstream file(filename, ifstream::in);
	string line;
	if (file.is_open()) {
		stack<mat4> modelview;
		modelview.push(mat4(1.0));
		vector<vec3> verts;
		vector<vec3> normverts;
		vector<vec3> norms;
		while (file.good()) {
			getline(file, line);
			parseLine(line, modelview, verts, normverts, norms);
		}
	} else {
		cerr << "Unable to open file " << filename << endl;
		exit(1);
	}
	file.close();
}
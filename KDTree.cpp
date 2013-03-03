/* Implementation of KDTree acceleration structure */

#include <iostream>
#include <vector>

#include "Shapes.h"
#include "Intersection.h"
#include "KDTree.h"

using namespace std;

#define TREELIMIT 5
#define EPSILON 0.000001

/***  TREENODE  ***/
struct ShapeSorter {
	ShapeSorter(int ax){this->axis = ax;}
	bool operator() (Shape* i, Shape* j) {
		return i->aabb.center[axis] < j->aabb.center[axis];
	}
	int axis;
};

TreeNode::TreeNode(vector<Shape*>& prims, int ax, AABB& bigBox, bool prevSize){
	this->aabb = bigBox;
	this->axis = ax;
	
	ShapeSorter s(axis);
	sort(prims.begin(), prims.end(),s);
	double split = prims[prims.size()/2]->aabb.center[axis];
	double leftSplit = DBL_MIN;
	double rightSplit = DBL_MAX;
	
	vector<Shape*> leftPrims;
	vector<Shape*>::iterator it=prims.begin();
	for(;it!=prims.end(); ++it){
		if((*it)->aabb.min[axis]<(split+EPSILON)){
			leftPrims.push_back(*it);
		}
	}
	
	vector<Shape*> rightPrims;
	for(it=prims.begin(); it!=prims.end(); ++it){
		if((*it)->aabb.max[axis]>(split-EPSILON)){
			rightPrims.push_back(*it);
		}
	}
	
	bool same = leftPrims.size()==prims.size() || rightPrims.size()==prims.size();
	
	//cout <<prims.size() << " "<< leftPrims.size() << " " << rightPrims.size() <<endl;
	
	if( prevSize&&same || prims.size()<=TREELIMIT){
		left = NULL;
		right = NULL;
		primatives = prims;
		//cout << "Leaf" << endl;
	}else{
		AABB leftAABB = AABB(aabb);
		leftAABB.max[axis] = split;
		
		AABB rightAABB = AABB(aabb);
		rightAABB.min[axis] = split;
		
		int newAxis = (axis+1)%3;
		left = new TreeNode(leftPrims, newAxis, leftAABB, same);
		right = new TreeNode(rightPrims, newAxis, rightAABB, same);
	}
	
}

TreeNode::~TreeNode(){
	if(left){
		delete left;
		delete right;
	}
}

Intersection TreeNode::intersect(Ray& ray){
	if(!left){
		return Intersection(primatives,ray);
	}
	
	double hitLeft = left->aabb.intersect(ray,this->axis);
	double hitRight = right->aabb.intersect(ray,this->axis);
	
	if(!hitLeft && !hitRight) return Intersection();
	
	if (!hitLeft) {
		return right->intersect(ray);
	}
	if (!hitRight) {
		return left->intersect(ray);
	}
	
	//hits both at same spot
	
	/*
	Intersection hit;
	if (abs(hitLeft-hitRight) < EPSILON){
		if (ray.direction[this->axis] > 0) {
			hit = right->intersect(ray);
			return hit;
		} else {
			hit = left->intersect(ray);
			return hit;
		}
	}
	*/
	
	Intersection hit;
	if (hitLeft > hitRight){
		hit = right->intersect(ray);
		if (!hit.primative){
			hit = left->intersect(ray);
		}
		return hit;
	} else {
		hit = left->intersect(ray);
		if (!hit.primative){
			hit = right->intersect(ray);
		}
		return hit;
	}
	
	/*
	
	Intersection hit1 = left->intersect(ray);
	Intersection hit2 = right->intersect(ray);
	if(!hit1.primative) return hit2;
	if(!hit2.primative) return hit1;
	if(glm::distance(hit1.point,ray.origin)<glm::distance(hit2.point,ray.origin)){
		return hit1;
	}
	return hit2;
	*/
	
	
}

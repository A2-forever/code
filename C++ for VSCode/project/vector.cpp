#include <cmath>
#include "vector.h"
using namespace std;

namespace VECTOR
{
	const double Rad_to_deg = 45.0 / atan(1.0);


	//private methods
	//set x by polar coordinate
	void Vector::set_x()
	{
		x = mag * cos(ang);
	}

	//set y by polar coordinate
	void Vector::set_y()
	{
		y = mag * sin(ang);
	}

	//set magnitude by Plane Rectangular Coordinate 
	void Vector::set_mag()
	{
		mag = sqrt(x * x + y * y);
	}

	//set angular by Plane Rectangular Coordinate 
	void Vector::set_ang()
	{
		if (x == 0 && y == 0)
			ang = 0.0;
		else
			ang = atan2(y, x);
	}


	//public methods
	//default constructor
	Vector::Vector()
	{
		x = y = ang = mag = 0;
		mode = RECT;
	}

	//construct vector from rectangular coordinates if form is RECT
	//or else form polar coordinates if form is POL
	Vector::Vector(double n1, double n2, Mode form)
	{
		mode = form;
		if (form == RECT)
		{
			x = n1;
			y = n2;
			set_mag();
			set_ang();
		}
		else if (form == POL)
		{
			mag = n1;
			ang = n2;
			set_x();
			set_y();
		}
		else
		{
			cout << "Incorrect 3 rd argument to Vector() -- ";
			cout << "Vector set to 0. " << endl;
			x = y = mag = ang = 0;
			mode = RECT;
		}
	}

	//reset vector from rectangular coordinates if form is RECT
	//or else form polar coordinates if form is POL
	 void Vector::reset(double n1, double n2, Mode form)
	{
		mode = form;
		if (form == RECT)
		{
			x = n1;
			y = n2;
			set_mag();
			set_ang();
		}
		else if (form == POL)
		{
			mag = n1;
			ang = n2;
			set_x();
			set_y();
		}
		else
		{
			cout << "Incorrect 3 rd argument to Vector() -- ";
			cout << "Vector set to 0. " << endl;
			x = y = mag = ang = 0;
			mode = RECT;
		}
	}

	 //destructor
	 Vector::~Vector()
	 {
	 }

	 //set to rectangular mode
	 void Vector::rect_mode()
	 {
		 mode = RECT;
	 }

	 //set to polar mode
	 void Vector::polar_mode()
	 {
		 mode = POL;
	 }


	 //operator overloading
	 //add two Vectors
	 Vector Vector::operator+(const Vector& b) const
	 {
		 return Vector(x + b.x, y + b.y);
	 }

	 //subsreact Vector b from a
	 Vector Vector::operator-(const Vector& b) const
	 {
		 return Vector(x - b.x, y - b.y);
	 }

	 //reverse sign of Vector
	 Vector Vector::operator-() const
	 {
		 return Vector(-x, -y);
	 }

	 //multiply vector by n
	 Vector Vector::operator*(double n) const
	 {
		 return Vector(n * x, n * y);
	 }


	 //friend methods
	 //multiply n by vector a
	 Vector operator* (double n, const Vector& a)
	 {
		 return Vector(n * a.x, n * a.y);
	 }

	 //display rectangular coordinates if mode is RECT,
	 //else display polar coordinates if mode is POL;
	 std::ostream& operator<<(std::ostream& os, const Vector& v)
	 {
		 if (v.mode == Vector::POL)
			 os << "(x,y)=(" << v.x << "," << v.y << ")";
		 else if (v.mode == Vector::RECT)
		 {
			 os << "(m,a)=(" << v.mag << "," << v.ang * Rad_to_deg << ")";
		 }
		 else
			 os << "Vector object mode is invalid";
		 return os;
	 }
}
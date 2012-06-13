/*
This library provides two quite distinct levels of usage. One is to use the basic class template interval<T> without specifying the policy. This only requires one to know and understand the concepts developed above and the content of the namespace boost. In addition to the class interval<T>, this level of usage provides arithmetic operators (+, -, *, /), algebraic and piecewise-algebraic functions (abs, square, sqrt, pow), transcendental and trigonometric functions (exp, log, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh), and the standard comparison operators (<, <=, >, >=, ==, !=), as well as several interval-specific functions (min, max, which have a different meaning than std::min and std::max; lower, upper, width, median, empty, singleton, equal, in, zero_in, subset, proper_subset, overlap, intersect, hull, bisect).
*/
#include <iostream>
#include <boost/numeric/interval.hpp>

int main(){
	boost::numeric::interval<float> t=0;
	for(int i=0;i<100000;i++)
		t+=0.1;
	std::cout<<lower(t)<<" "<<upper(t)<<" "<<width(t)<<std::endl;
}

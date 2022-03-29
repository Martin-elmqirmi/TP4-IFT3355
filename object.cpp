#include "object.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>


bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assure une valeur correcte pour la coordonnée W de l'origine et de la direction
	// Vous pouvez commentez ces lignes si vous faites très attention à la façon de construire vos rayons.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray local_ray(i_transform * ray.origin, i_transform * ray.direction);
	//!!! NOTE UTILE : pour calculer la profondeur dans localIntersect(), si l'intersection se passe à
	// ray.origin + ray.direction * t, alors t est la profondeur
	//!!! NOTE UTILE : ici, la direction peut êytre mise à l'échelle, alors vous devez la renormaliser
	// dans localIntersect(), ou vous aurez une profondeur dans le système de coordonnées local, qui
	// ne pourra pas être comparée aux intersection avec les autres objets.
    if (localIntersect(local_ray, hit)) 
	{
        // Assure la valeur correcte de W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transforme les coordonnées de l'intersection dans le repère global.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}


bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
    // @@@@@@ VOTRE CODE ICI
	// Vous pourriez aussi utiliser des relations géométriques pures plutôt que les
	// outils analytiques présentés dans les slides.
	// Ici, dans le système de coordonées local, la sphère est centrée en (0, 0, 0)
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
    Vector o = ray.origin;
    Vector d = ray.direction;
    double A = d.dot(d);
    double B = 2 * (o.dot(d));
    double C = o.dot(o) - pow(radius, 2);
    double D = pow(B, 2) - 4 * A * C;

    if(D < 0) {
        return false;
    } else if(D > -0.000001f and D < 0.0000001f) {
        double t = - B / (2 * A);
        // L'intersection est à l'inverse de la direction du rayon
        if(t > 0) {
            if(hit.depth > t) {
                hit.depth = t;
                hit.position = o + t * d;
                hit.normal = - hit.position.normalized();
                return true;
            } else return false;
        } else return false;
    } else {
        double t1 = - (B + sqrt(D)) / (2 * A);
        double t2 = - (B - sqrt(D)) / (2 * A);
        if(t1 < 0 and t2 < 0) { // les intersections sont à l'inverse de la direction du rayon
            return false;
        } else if(t1 < 0 or t2 < 0) { // Le rayon est à l'intérieur de la sphère
            double tmax = std::max(t1, t2);
            if(hit.depth > tmax) {
                hit.depth = tmax;
                hit.position = o + tmax * d;
                hit.normal = - hit.position.normalized();
            } else return false;
        } else { // Le rayon intersecte la sphère à 2 endroits
            // Je choisis la première intersection
            double tmin = std::min(t1, t2);
            if(hit.depth > tmin) {
                hit.depth = tmin;
                hit.position = o + tmin * d;
                hit.normal = - hit.position.normalized();
            } else return false;
        }
        return true;
    }
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
	// @@@@@@ VOTRE CODE ICI
	// N'acceptez pas les intersections tant que le rayon est à l'intérieur du plan.
	// ici, dans le système de coordonées local, le plan est à z = 0.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
    Vector o = ray.origin;
    Vector d = ray.direction;
    Vector normalPlan = Vector(0.0, 0.0, - 1.0);
    double dot = d.dot(normalPlan);
    if(std::abs(dot) < 0.00001f) {
        return false;
    }
    double t = - o[2] / d[2];
    if(t < 0 or hit.depth < t) {
        return false;
    }

    hit.depth = t;
    hit.position = o + t * d;
    hit.normal = normalPlan.normalized();

    return true;
}


bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
    // @@@@@@ VOTRE CODE ICI (licence créative)
    return false;
}


// Intersections !
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Test de la boite englobante
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Rayon parallèle à un plan de la boite englobante et en dehors de la boite
				return false;
			}
			// Rayon parallèle à un plan de la boite et dans la boite: on continue
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Assure t1 <= t2

			if (t1 > tNear) tNear = t1; // On veut le plus lointain tNear.
			if (t2 < tFar) tFar = t2; // On veut le plus proche tFar.

			if (tNear > tFar) return false; // Le rayon rate la boite englobante.
			if (tFar < 0) return false; // La boite englobante est derrière le rayon.
		}
	}
	// Si on arrive jusqu'ici, c'est que le rayon a intersecté la boite englobante.

	// Le rayon interesecte la boite englobante, donc on teste chaque triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const &tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}

bool Mesh::intersectTriangle(Ray const &ray,
	Triangle const &tri,
	Intersection &hit) const
{
	// Extrait chaque position de sommet des données du maillage.
	Vector const &p0 = positions[tri[0].pi];
	Vector const &p1 = positions[tri[1].pi];
	Vector const &p2 = positions[tri[2].pi];

	// @@@@@@ VOTRE CODE ICI
	// Décidez si le rayon intersecte le triangle (p0,p1,p2).
	// Si c'est le cas, remplissez la structure hit avec les informations
	// de l'intersection et renvoyez true.
	// Vous pourriez trouver utile d'utiliser la routine implicitLineEquation()
	// pour calculer le résultat de l'équation de ligne implicite en 2D.
	//
	// NOTE : hit.depth est la profondeur de l'intersection actuellement la plus proche,
	// donc n'acceptez pas les intersections qui occurent plus loin que cette valeur.
	//!!! NOTE UTILE : pour le point d'intersection, sa normale doit satisfaire hit.normal.dot(ray.direction) < 0
    Vector normal = ((p2 - p0).cross(p1 - p0)).normalized();
    Vector p = ray.origin;
    Vector r = ray.direction;

    double D = - (p0.dot(normal));
    double t = (- normal.dot(p) - D) / (normal.dot(r));

    if(t <= 0 or hit.depth < t) {
        return false;
    }

    Vector v01 = p1 - p0;
    Vector v12 = p2 - p1;
    Vector v20 = p0 - p2;

    Vector position = p + t * r;
    Vector v0p = position - p0;
    Vector v1p = position - p1;
    Vector v2p = position - p2;

    Vector alpha = v0p.cross(v01);
    Vector beta = v1p.cross(v12);
    Vector gamma = v2p.cross(v20);

    double dotAB = alpha.dot(beta);
    double dotAG = alpha.dot(gamma);
    double dotBG = beta.dot(gamma);

    if(!((dotAB >= 0 and dotAG >= 0 and dotBG >= 0) or (dotAB < 0 and dotAG < 0 and dotBG < 0))){
        return false;
    }

    hit.depth = t;
    hit.position = position;
    hit.normal = normal;

    return true;
}

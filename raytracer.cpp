#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

#include "raytracer.hpp"
#include "image.hpp"


void Raytracer::render(const char *filename, const char *depth_filename,
                       Scene const &scene)
{
    // Alloue les deux images qui seront sauvegardées à la fin du programme.
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);
    
    // Crée le zBuffer.
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        zBuffer[i] = DBL_MAX;
    }

	// @@@@@@ VOTRE CODE ICI
	// Calculez les paramètres de la caméra pour les rayons. Référez-vous aux slides pour les détails.
	//!!! NOTE UTILE : tan() prend des radians plutot que des degrés. Utilisez deg2rad() pour la conversion.
	//!!! NOTE UTILE : Le plan de vue peut être n'importe où, mais il sera implémenté différement.
	// Vous trouverez des références dans le cours.
    Vector v = - scene.camera.up.normalized();
    Vector w = (scene.camera.center - scene.camera.position).normalized();
    Vector u = v.cross(w).normalized();
    Vector C = scene.camera.position;
    double d = scene.camera.zNear;
    double t = tan(deg2rad(scene.camera.fovy / 2)) * d;
    double b = - t;
    double r = t * scene.camera.aspect;
    double l = - r;

    // Itère sur tous les pixels de l'image.
    for(int y = 0; y < scene.resolution[1]; y++) {
        for(int x = 0; x < scene.resolution[0]; x++) {

            // Génère le rayon approprié pour ce pixel.
			Ray ray;
			if (scene.objects.empty())
			{
				// Pas d'objet dans la scène --> on rend la scène par défaut.
				// Pour celle-ci, le plan de vue est à z = 640 avec une largeur et une hauteur toute deux à 640 pixels.
				ray = Ray(scene.camera.position, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - scene.camera.position).normalized());
			}
			else
			{
				// @@@@@@ VOTRE CODE ICI
				// Mettez en place le rayon primaire en utilisant les paramètres de la caméra.
				//!!! NOTE UTILE : tous les rayons dont les coordonnées sont exprimées dans le
				//                 repère monde doivent avoir une direction normalisée.
                Vector O = C + (d * w) + (l * u) + (t * v);
                Vector Pij = O + ((x + 0.5) * (r - l) / scene.resolution[0]) * u - ((y + 0.5) * (t - b) / scene.resolution[1]) * v;
                ray = Ray(C, (Pij - C).normalized());
			}

            // Initialise la profondeur de récursivité du rayon.
            int rayDepth = 0;
           
            // Notre lancer de rayons récursif calculera la couleur et la z-profondeur.
            Vector color;

            // Ceci devrait être la profondeur maximum, correspondant à l'arrière plan.
            // NOTE : Ceci suppose que la direction du rayon est de longueur unitaire (normalisée)
			//        et que l'origine du rayon est à la position de la caméra.
            double depth = scene.camera.zFar;

            // Calcule la valeur du pixel en lançant le rayon dans la scène.
            trace(ray, rayDepth, scene, color, depth);

            // Test de profondeur
            if(depth >= scene.camera.zNear && depth <= scene.camera.zFar &&
                depth < zBuffer[x + y*scene.resolution[0]]) {
                zBuffer[x + y*scene.resolution[0]] = depth;

                // Met à jour la couleur de l'image (et sa profondeur)
                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) /
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

		// Affiche les informations de l'étape
		if (y % 100 == 0)
		{
			printf("Row %d pixels finished.\n", y);
		}
    }

	// Sauvegarde l'image
    colorImage.writeBMP(filename);
    depthImage.writeBMP(depth_filename);

	printf("Ray tracing finished with images saved.\n");

    delete[] zBuffer;
}


bool Raytracer::trace(Ray const &ray, 
                 int &rayDepth,
                 Scene const &scene,
                 Vector &outColor, double &depth)
{
    // Incrémente la profondeur du rayon.
    rayDepth++;

    // - itérer sur tous les objets en appelant   Object::intersect.
    // - ne pas accepter les intersections plus lointaines que la profondeur donnée.
    // - appeler Raytracer::shade avec l'intersection la plus proche.
    // - renvoyer true ssi le rayon intersecte un objet.
	if (scene.objects.empty())
	{
		// Pas d'objet dans la scène --> on rend la scène par défaut :
		// Par défaut, un cube est centré en (0, 0, 1280 + 160) avec une longueur de côté de 320, juste en face de la caméra.
		// Test d'intersection :
		double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
		double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
		if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
		{
			// S'il y a intersection :
			Material m; m.emission = Vector(16.0, 0, 0); m.reflect = 0; // seulement pour le matériau par défaut ; vous devrez utiliser le matériau de l'objet intersecté
			Intersection intersection;	// seulement par défaut ; vous devrez passer l'intersection trouvée par l'appel à Object::intersect()
			outColor = shade(ray, rayDepth, intersection, m, scene);
			depth = 1280;	// la profondeur devrait être mise à jour dans la méthode Object::intersect()
		}
	}
	else
	{
		// @@@@@@ VOTRE CODE ICI
		// Notez que pour Object::intersect(), le paramètre hit correspond à celui courant.
		// Votre intersect() devrait être implémenté pour exclure toute intersection plus lointaine que hit.depth
        Intersection intersection;
        intersection.depth = scene.camera.zFar;
        Material m;
        bool intersect = false;
        for (auto& object : scene.objects) {
            if(object->intersect(ray, intersection)) {
                m = object->material;
                outColor = shade(ray, rayDepth, intersection, m, scene);
                depth = intersection.depth;
                intersect = true;
            }
        }
        return intersect;
	}

    // Décrémente la profondeur du rayon.
    rayDepth--;

    return false; 
}


Vector Raytracer::shade(Ray const &ray,
                 int &rayDepth,
                 Intersection const &intersection,
                 Material const &material,
                 Scene const &scene)
{
    // - itérer sur toutes les sources de lumières, calculant les contributions ambiant/diffuse/speculaire
    // - utiliser les rayons d'ombre pour déterminer les ombres
    // - intégrer la contribution de chaque lumière
    // - inclure l'émission du matériau de la surface, s'il y a lieu
    // - appeler Raytracer::trace pour les couleurs de reflection/refraction
    // Ne pas réfléchir/réfracter si la profondeur de récursion maximum du rayon a été atteinte !
	//!!! NOTE UTILE : facteur d'atténuation = 1.0 / (a0 + a1 * d + a2 * d * d)..., la lumière ambiante ne s'atténue pas, ni n'est affectée par les ombres
	//!!! NOTE UTILE : n'acceptez pas les intersection des rayons d'ombre qui sont plus loin que la position de la lumière
	//!!! NOTE UTILE : pour chaque type de rayon, i.e. rayon d'ombre, rayon reflechi, et rayon primaire, les profondeurs maximales sont différentes
	Vector diffuse(0);
	Vector ambient(0);
	Vector specular(0);
	for (auto lightIter = scene.lights.begin(); lightIter != scene.lights.end(); lightIter++)
	{
		// @@@@@@ VOTRE CODE ICI
		// Calculez l'illumination locale ici, souvenez-vous d'ajouter les lumières ensemble.
		// Testez également les ombres ici, si un point est dans l'ombre aucune couleur (sauf le terme ambient) ne devrait être émise.
        // Direction de la lumiere
        Vector lightDirection = - (intersection.position - lightIter->position).normalized();
        Vector viewDirection = (intersection.position - scene.camera.position).normalized();
        Vector intersectionNormal = intersection.normal.normalized();

        // Ambient
        ambient += lightIter->ambient * material.ambient;

        // Shadows
        Vector lightRayDirection = lightIter->position - intersection.position;
        double lightDistance = lightRayDirection.length();
        Vector positionAcnee = intersection.position + 0.001 * lightRayDirection;
        Ray lightRay = Ray(positionAcnee, lightRayDirection.normalized());

        Intersection inter;
        inter.depth = scene.camera.zFar;
        bool notIntersect = true;
        for (auto& object : scene.objects) {
            if(object->intersect(lightRay, inter)){
                notIntersect = false;
            };
        }

        // Si on intersecte un objet après la lumière on est pas dans l'ombre
        if(notIntersect or lightDistance < inter.depth) {
            // Facteur d'atténuation;
            Vector att = lightIter->attenuation;
            double attenuation = 1.0 / (att[0] + att[1] * lightDistance + att[2] * pow(lightDistance, 2));

            // Diffuse
            double diffuseDot = std::max(intersectionNormal.dot(lightDirection), 0.0);
            diffuse += lightIter->diffuse * material.diffuse * diffuseDot * attenuation;

            // Specular
            Vector h = ((lightDirection + viewDirection) / 2.0).normalized();
            double specularDot = std::max((h).dot(intersectionNormal), 0.0);
            specular += lightIter->specular * material.specular * pow(specularDot, material.shininess) * attenuation;
        }
	}

	Vector reflectedLight(0);
	if ((!(ABS_FLOAT(material.reflect) < 1e-6)) && (rayDepth < 100))
	{
		// @@@@@@ VOTRE CODE ICI
		// Calculez la couleur réfléchie en utilisant trace() de manière récursive.
        Vector color;
        double depth = scene.camera.zFar;
        Vector rayDirection = ray.direction.normalized();
		Vector newRayDirection = (rayDirection + 2 * (intersection.normal.dot(- rayDirection)) * intersection.normal).normalized();
        Ray newRay = Ray(intersection.position + 0.000001 * newRayDirection, newRayDirection);
        trace(newRay, rayDepth, scene, color, depth);
        reflectedLight += color;
	}

	return material.emission + ambient + diffuse + specular + material.reflect * reflectedLight;
}
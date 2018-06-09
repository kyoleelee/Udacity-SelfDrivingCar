/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

default_random_engine gen;

static int NUM_PARTICLES = 200;

//for the update step
double bivariate_normal(double x, double y, double mu_x, double mu_y, double sig_x, double sig_y) {
	return exp(-((x-mu_x)*(x-mu_x)/(2*sig_x*sig_x) + (y-mu_y)*(y-mu_y)/(2*sig_y*sig_y))) / (2.0*3.14159*sig_x*sig_y);
}



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Set number of particles and weight size
	num_particles = NUM_PARTICLES;
	particles.resize(num_particles);
	weights.resize(num_particles);
	
	// Create Gaussin distribution noise
	normal_distribution<double> x_dist(x, std[0]);
	normal_distribution<double> y_dist(y, std[1]);
	normal_distribution<double> theta_dist(theta, std[2]);

	// Init particles
	for (int i=0; i<num_particles; i++) {
		Particle tmp;
		
		tmp.id = i;
		tmp.x = x_dist(gen);
		tmp.y = y_dist(gen);
		tmp.theta = theta_dist(gen);
		tmp.weight = 1;

		particles.push_back(tmp);
	}
		
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	for (int i=0; i<num_particles; i++) {
		// Generate Gussian distribution noise
		normal_distribution<double> x_dist(0, std_pos[0]);
		normal_distribution<double> y_dist(0, std_pos[1]);
		normal_distribution<double> theta_dist(0, std_pos[2]);
		
		if (fabs(yaw_rate) < 0.0001) {			
			// Yaw is constant
			// Update the particle position 			
			particles[i].x += velocity*cos(particles[i].theta)*delta_t;
			particles[i].y += velocity*sin(particles[i].theta)*delta_t;
		}

		else {
			// Yaw is not constant			
			// Update the particle position
			particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;
		}

		// Set particle state
		particles[i].x += x_dist(gen);
		particles[i].y += y_dist(gen);
		particles[i].theta += theta_dist(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for (unsigned int i=0; i<observations.size(); i++) {
		double min_dist = numeric_limits<float>::max();

		for (unsigned int j=0; j<predicted.size(); j++) {
			double distance = dist(observations[i].x,observations[i].y, predicted[j].x, predicted[j].y);

			if (distance < min_dist) {
				min_dist = distance;
				observations[i].id = predicted[j].id;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	// Clear weights
	weights.clear();


  for (unsigned int i=0; i<particles.size(); i++) {
    particles[i].weight = 1.0;

    // step 1: collect valid landmarks
    vector<LandmarkObs> predictions;

    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); j++) {
      double distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);

      if (distance < sensor_range){
				// if the landmark is within the sensor range, save it to predictions
      	LandmarkObs pred_landmark;
				pred_landmark.id = map_landmarks.landmark_list[j].id_i;
				pred_landmark.x = map_landmarks.landmark_list[j].x_f;  
				pred_landmark.y = map_landmarks.landmark_list[j].y_f;

				predictions.push_back(pred_landmark);
      }
    }

    // step 2: convert observations coordinates from vehicle to map
    vector<LandmarkObs> observations_map;

		for (unsigned int j=0; j<observations.size(); j++){
			LandmarkObs m_observation;

			m_observation.x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
			m_observation.y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
			m_observation.id = -1;

			observations_map.push_back(m_observation);
		}

    // step 3: find landmark index for each observation
    dataAssociation(predictions, observations_map);

    // step 4: compute the particle's weight:
    for (unsigned int l=0; l<observations_map.size(); l++) {

      Map::single_landmark_s landmark = map_landmarks.landmark_list.at(observations_map[l].id-1);
    	// see equation this link:
      double x_term = pow(observations_map[l].x - landmark.x_f, 2) / (2 * pow(std_landmark[0], 2));
      double y_term = pow(observations_map[l].y - landmark.y_f, 2) / (2 * pow(std_landmark[1], 2));
      double w = exp(-(x_term + y_term)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      particles[i].weight *= w;
    }

    weights.push_back(particles[i].weight);
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	random_device rd_wts;
	mt19937 generator_wts(rd_wts());

	// Creates a discrete distribution for weight.
	discrete_distribution<int> distribution_wts(weights.begin(), weights.end());
	vector<Particle> resampled_particles;

	// Resample
	for (int i=0; i<num_particles; i++) {
		Particle tmp = particles[distribution_wts(generator_wts)];
		resampled_particles.push_back(tmp);
	}

	particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

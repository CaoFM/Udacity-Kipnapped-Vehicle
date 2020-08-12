/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  num_particles = 1000;  // TODO: Set the number of particles
  
  for (int i = 0; i<num_particles; ++i) {
    Particle sample;
   
    sample.id = i;
    sample.x = dist_x(gen);
    sample.y = dist_y(gen);
    sample.theta = dist_theta(gen);
    sample.weight = 1.0;

    particles.push_back(sample);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;
  
  for (int i=0;i<particles.size();++i) {
    if (yaw_rate == 0) {
      particles[i].x += (velocity * sin(particles[i].theta)*delta_t);
      particles[i].y += (velocity * cos(particles[i].theta)*delta_t);
    }
    else {
      particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    /* normalize theta to 0 ~ 2pi */
    while (particles[i].theta >= (2*M_PI)){
      particles[i].theta -= (2*M_PI);
    }
    while (particles[i].theta < 0){
      particles[i].theta += (2*M_PI);
    }
 
  }
  /* add Gaussian to each particle */  

  for (int i=0;i<particles.size();++i) {
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  /* transform vehicle coordinate observation to map coordinate */
  for (int p=0; p<particles.size(); ++p) {
    /* for each particle calculate its weight */
    double temp_weight = 1.0;
    particles[p].associations.clear();
    particles[p].sense_x.clear();    
    particles[p].sense_y.clear();    
    for (int o=0; o<observations.size(); ++o) {
      /* for each observation translate to map coordinate */
      double x_m = particles[p].x 
                   + (cos(particles[p].theta) * observations[o].x)
                   - (sin(particles[p].theta) * observations[o].y);
      double y_m = particles[p].y 
                   + (sin(particles[p].theta) * observations[o].x)
                   + (cos(particles[p].theta) * observations[o].y);
      
      /* find out the closest land mark to world coordinate [x_m, y_m] */
      int nearest_id;
      double distance;
      double min_distance;
      double landmark_x;
      double landmark_y;
      
      for (int m=0; m<map_landmarks.landmark_list.size(); ++m) {
        /* calculate distance between observation [o] and landmark [m]*/
        distance = dist(x_m,
                        y_m,
                        map_landmarks.landmark_list[m].x_f,
                        map_landmarks.landmark_list[m].y_f);
        
        if (   (m == 0)
             ||(   (m>0)
                && (distance < min_distance)
              )
          )
        {
          /* first map landmark - initialize the closest fit */
          /* afterwards         - update if distance is shorter */
          landmark_x = map_landmarks.landmark_list[m].x_f;
          landmark_y = map_landmarks.landmark_list[m].y_f;
          nearest_id = map_landmarks.landmark_list[m].id_i;
          min_distance = distance;
        }
      }
      
      /* one observation is now tranformed and associated to a landmark id -- register it */
      particles[p].associations.push_back(nearest_id);
      particles[p].sense_x.push_back(x_m);    
      particles[p].sense_y.push_back(y_m);
 
      /* calculate weight of the particle to the associated landmark */
      /* the associated landmark has been recorded during association
         and have the x, y value stored in temp variable of landmark_x and landmark_y */
      double prob =( ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])
                      *exp(-(   (pow((landmark_x-x_m),2)/(2*std_landmark[0]*std_landmark[0]))
                               +(pow((landmark_y-y_m),2)/(2*std_landmark[1]*std_landmark[1]))
                             )
                           )
                      ));
      
      temp_weight *= prob;
    }
    /* temporarily store the raw weight */
    particles[p].weight = temp_weight;
  } 

  /* sum weight from all particles */
  double sum_weight = 0;
  for (int p=0; p<particles.size(); ++p) {
    sum_weight += particles[p].weight;
  }
  
  /* normalize weight */
  for (int p=0; p<particles.size(); ++p) {
    particles[p].weight /= sum_weight;
  }

}

void ParticleFilter::resample() {
   /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::vector<double> w;
  
  for (int i = 0; i < particles.size(); ++i) {
    w.push_back(particles[i].weight);
  }

  std::discrete_distribution<int> distribution(w.begin(),w.end());
  
  std::vector<Particle> resampled;
  
  for (int i=0; i<particles.size(); ++i) {
    int number = distribution(gen);
    resampled.push_back(particles[number]);
  }
  
  particles = resampled;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
#include "hermes2d.h"

#define PID_DEFAULT_TOLERANCE 0.25
#define DEFAULT_STEP 0.1
#define STEP_TOLERANCE 1.0e-10

//HERMES_API_USED_TEMPLATE(Hermes::vector<Solution*>);

class HERMES_API PidTimestepController {

public:
  PidTimestepController(double final_time, bool pid_on = true,
      double default_step = DEFAULT_STEP, double scaling=1.0, double tolerance = PID_DEFAULT_TOLERANCE, bool err_out=true) {
    this->delta = tolerance;
    this->final_time = final_time/scaling;
    this->time = 0;
    this->step_number = 0;
    this->scaling = scaling;
    timestep = new double;
    (*timestep) = default_step/scaling;
    this->pid = pid_on;
    this->err_out = err_out;
    finished = false;
  };

  ~PidTimestepController() {
    if (graphs != NULL) {
      //delete[] graphs;
    }
  }

  int get_timestep_number() {return step_number;};
  double get_time() {return time*scaling;};
  double get_scaled_time() {return time;};

  // true if next time step can be run, false if the time step must be re-run with smaller time step.
  bool end_step(Hermes::vector<Solution*> solutions, Hermes::vector<Solution *> prev_solutions);
  void begin_step();
  bool has_next();

  // reference to the current calculated time step
  double *timestep;

private:
  bool pid;
  bool err_out;
  double delta;
  double scaling;
  double final_time;
  double time;
  int step_number;
  std::vector<double> err_vector;
  bool finished;
  SimpleGraph* graphs;

  // PID parameters
  const static double kp = 0.075;
  const static double kl = 0.175;
  const static double kD = 0.01;

};


// Usage: do{ begin_step() calculations .... end_step(..);} while(has_next());

void PidTimestepController::begin_step() {

  if ((time + (*timestep)) >= final_time
      || std::abs(time + (*timestep) - final_time) < STEP_TOLERANCE) {
    info("Time step would exceed the final time... reducing");
    (*timestep) = final_time - time;
    info("The last time step: %g", *timestep);
    finished = true;
  }
  time += (*timestep);
  step_number++;
  info("begin_step processed, new step number: %i and cumulative time: %g", step_number, time);
}

bool PidTimestepController::end_step(Hermes::vector<Solution *> solutions,
    Hermes::vector<Solution *> prev_solutions) {

  if (pid) {

    unsigned int neq = solutions.size();
    if (neq == 0) {
      return true;
    }
    if (prev_solutions.size() == 0) {
      return true;
    }
    if (neq != prev_solutions.size()) {
      error_function("Inconsistent parameters in PidTimestepController::next(...)");
    }
    double max_rel_error = 0.0;

    for (unsigned int i = 0; i < neq; i++) {
      double rel_error = calc_rel_error(solutions[i], prev_solutions[i], HERMES_H1_NORM);
      max_rel_error = (rel_error > max_rel_error) ? rel_error : max_rel_error;
      if (err_out) {
        if (graphs == NULL) {
          graphs = new SimpleGraph[neq];
        }

        graphs[i].add_values(get_time()*scaling, rel_error);
        graphs[i].save_numbered("time_relerr%i.dat", (int)i);

      }
      info("Solution[%i]: rel error %g, largest relative error %g",
          i, rel_error, max_rel_error);
    }

    err_vector.push_back(max_rel_error);

    if (err_vector.size() > 2 && max_rel_error <= delta) {
      int size = err_vector.size();
      info("Error vector sufficient for adapting...");
      double t_coeff = pow(err_vector.at(size - 2)/err_vector.at(size-1),kp)
          * pow(0.25/err_vector.at(size - 1), kl)
          * pow(err_vector.at(size - 2)*err_vector.at(size - 2)/(err_vector.at(size -1)*err_vector.at(size-3)), kD);
       info("Coefficient %g", t_coeff);
       (*timestep) = (*timestep)*t_coeff;
       info("New time step: %g", *timestep);
    }
  } // end pid

  return true;

}

bool PidTimestepController::has_next() {
  return !finished;
}



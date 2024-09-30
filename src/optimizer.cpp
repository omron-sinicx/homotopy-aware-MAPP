#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_multi_edge.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/sparse_optimizer_terminate_action.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>

#include <yaml-cpp/yaml.h>

const double EPS = 1e-9;

class WayPointVertex : public g2o::BaseVertex<2, Eigen::Vector2d> {
public:
  void setToOriginImpl() { _estimate.setZero(); }
  void oplusImpl(const double *update) {
    _estimate[0] += update[0];
    _estimate[1] += update[1];
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class ParameterVertex : public g2o::BaseVertex<1, double> {
public:
  ParameterVertex() {}
  void setToOriginImpl() { _estimate = 0.; }
  void oplusImpl(const double *update) { _estimate += *update; }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class CollisionConstraint
    : public g2o::BaseBinaryEdge<1, double, WayPointVertex, WayPointVertex> {
public:
  CollisionConstraint(WayPointVertex *point_0, WayPointVertex *point_1,
                      const double agent_radius, const double coefficient)
      : g2o::BaseBinaryEdge<1, double, WayPointVertex, WayPointVertex>() {
    _vertices[0] = point_0;
    _vertices[1] = point_1;
    _measurement = 2 * agent_radius;
    _information[0] = coefficient;
  }
  void setCoefficient(const double coefficient) {
    _information[0] = coefficient;
  }

  void computeError() {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    double distance = (point_1 - point_0).norm();
    if (distance >= _measurement) {
      _error[0] = 0.;
    } else {
      _error[0] = 1. / (distance + EPS) - 1. / _measurement;
    }
    assert(std::isfinite(_error[0]));
  }
  void linearizeOplus() override {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    double distance = (point_1 - point_0).norm();
    if (distance >= 1.0) {
      std::get<0>(_jacobianOplus) = Eigen::MatrixXd::Zero(1, 2);
      std::get<1>(_jacobianOplus) = Eigen::MatrixXd::Zero(1, 2);
    } else {
      double dedd = -std::pow(distance + EPS, -2);
      Eigen::Vector2d direction = (point_1 - point_0) / (distance + EPS);
      std::get<0>(_jacobianOplus) = (-dedd * direction).transpose();
      std::get<1>(_jacobianOplus) = (dedd * direction).transpose();
    }
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class PathCostConstraint : public g2o::BaseMultiEdge<2, Eigen::Vector2d> {
public:
  PathCostConstraint(WayPointVertex *point_0, WayPointVertex *point_1,
                     ParameterVertex *delta, const double coefficient)
      : g2o::BaseMultiEdge<2, Eigen::Vector2d>() {
    resize(3);
    _vertices[0] = point_0;
    _vertices[1] = point_1;
    _vertices[2] = delta;
    _information = coefficient * Eigen::Matrix2d::Identity();
  }
  void computeError() {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[2])->estimate();
    _error = (point_1 - point_0) / delta;
  }
  void linearizeOplus() override {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[2])->estimate();
    _jacobianOplus[0] = -1. / delta * Eigen::Matrix2d::Identity();
    _jacobianOplus[1] = 1. / delta * Eigen::Matrix2d::Identity();
    _jacobianOplus[2] = -std::pow(delta, -2) * (point_1 - point_0);
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class PathCostConstraintMaxVelocity : public g2o::BaseMultiEdge<1, double> {
public:
  PathCostConstraintMaxVelocity(WayPointVertex *point_0,
                                WayPointVertex *point_1, ParameterVertex *delta,
                                const double max_velocity,
                                const double coefficient)
      : g2o::BaseMultiEdge<1, double>() {
    resize(3);
    _vertices[0] = point_0;
    _vertices[1] = point_1;
    _vertices[2] = delta;
    _measurement = max_velocity;
    _information[0] = coefficient;
  }
  void setCoefficient(const double coefficient) {
    _information[0] = coefficient;
  }
  void computeError() {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[2])->estimate();
    _error[0] = std::max(0., (point_1 - point_0).norm() / delta - _measurement);
  }
  void linearizeOplus() override {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    double delta_ =
        static_cast<const ParameterVertex *>(_vertices[2])->estimate();
    double distance = (point_1 - point_0).norm();
    if (distance / delta_ < _measurement) {
      _jacobianOplus[0] = Eigen::MatrixXd::Zero(1, 2);
      _jacobianOplus[1] = Eigen::MatrixXd::Zero(1, 2);
      _jacobianOplus[2] = Eigen::MatrixXd::Zero(1, 1);
    } else {
      Eigen::Vector2d direction = (point_1 - point_0) / (distance + EPS);
      // double r = distance/delta - _measurement;
      _jacobianOplus[0] = (-direction / delta_).transpose();
      _jacobianOplus[1] = (direction / delta_).transpose();
      _jacobianOplus[2](0, 0) = -std::pow(delta_, -2) * distance;
    }
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class PathCostConstraintMaxAcceleration : public g2o::BaseMultiEdge<1, double> {
public:
  PathCostConstraintMaxAcceleration(WayPointVertex *point_0,
                                    WayPointVertex *point_1,
                                    WayPointVertex *point_2,
                                    ParameterVertex *delta,
                                    const double max_acceleration,
                                    const double coefficient)
      : g2o::BaseMultiEdge<1, double>() {
    resize(4);
    _vertices[0] = point_0;
    _vertices[1] = point_1;
    _vertices[2] = point_2;
    _vertices[3] = delta;
    _measurement = max_acceleration;
    _information[0] = coefficient;
  }
  void setCoefficient(const double coefficient) {
    _information[0] = coefficient;
  }
  void computeError() {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    Eigen::Vector2d point_2 =
        static_cast<const WayPointVertex *>(_vertices[2])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[3])->estimate();
    double acceleration =
        (point_0 + point_2 - 2. * point_1).norm() / std::pow(delta, 2);
    _error[0] = std::max(0., acceleration - _measurement);
  }
  void linearizeOplus() override {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    Eigen::Vector2d point_2 =
        static_cast<const WayPointVertex *>(_vertices[2])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[3])->estimate();
    double acceleration = (point_0 + point_2 - 2. * point_1).norm();
    if (acceleration / std::pow(delta, 2) < _measurement) {
      _jacobianOplus[0] = Eigen::MatrixXd::Zero(1, 2);
      _jacobianOplus[1] = Eigen::MatrixXd::Zero(1, 2);
      _jacobianOplus[2] = Eigen::MatrixXd::Zero(1, 2);
      _jacobianOplus[3] = Eigen::MatrixXd::Zero(1, 1);
    } else {
      Eigen::Vector2d direction =
          (point_0 + point_2 - 2. * point_1) / acceleration;
      _jacobianOplus[0] = std::pow(delta, -2) * direction.transpose();
      _jacobianOplus[1] = -2. * std::pow(delta, -2) * direction.transpose();
      _jacobianOplus[2] = std::pow(delta, -2) * direction.transpose();
      _jacobianOplus[3](0, 0) = -2. * std::pow(delta, -3) * acceleration;
    }
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class PathCostConstraintAcceleration
    : public g2o::BaseMultiEdge<2, Eigen::Vector2d> {
public:
  PathCostConstraintAcceleration(WayPointVertex *point_0,
                                 WayPointVertex *point_1,
                                 WayPointVertex *point_2,
                                 ParameterVertex *delta,
                                 const double coefficient)
      : g2o::BaseMultiEdge<2, Eigen::Vector2d>() {
    resize(4);
    _vertices[0] = point_0;
    _vertices[1] = point_1;
    _vertices[2] = point_2;
    _vertices[3] = delta;
    _information = coefficient * Eigen::Matrix2d::Identity();
  }
  void computeError() {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    Eigen::Vector2d point_2 =
        static_cast<const WayPointVertex *>(_vertices[2])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[3])->estimate();
    _error = (point_0 + point_2 - 2. * point_1) / std::pow(delta, 2);
  }
  void linearizeOplus() override {
    Eigen::Vector2d point_0 =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    Eigen::Vector2d point_1 =
        static_cast<const WayPointVertex *>(_vertices[1])->estimate();
    Eigen::Vector2d point_2 =
        static_cast<const WayPointVertex *>(_vertices[2])->estimate();
    double delta =
        static_cast<const ParameterVertex *>(_vertices[3])->estimate();
    _jacobianOplus[0] = std::pow(delta, -2) * Eigen::Matrix2d::Identity();
    _jacobianOplus[1] = -2. * std::pow(delta, -2) * Eigen::Matrix2d::Identity();
    _jacobianOplus[2] = std::pow(delta, -2) * Eigen::Matrix2d::Identity();
    _jacobianOplus[3] =
        -2. * std::pow(delta, -3) * (point_0 + point_2 - 2. * point_1);
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class ParameterTargetConstraint
    : public g2o::BaseUnaryEdge<1, double, ParameterVertex> {
public:
  ParameterTargetConstraint(ParameterVertex *parameter, const double target,
                            const double coefficient)
      : g2o::BaseUnaryEdge<1, double, ParameterVertex>() {
    _vertices[0] = parameter;
    _measurement = target;
    _information[0] = coefficient;
  }
  void computeError() {
    _error[0] = static_cast<const ParameterVertex *>(_vertices[0])->estimate() -
                _measurement;
  }
  void linearizeOplus() override { std::get<0>(_jacobianOplus)(0, 0) = 1.; }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class WayPointTargetConstraint
    : public g2o::BaseUnaryEdge<2, Eigen::Vector2d, WayPointVertex> {
public:
  WayPointTargetConstraint(WayPointVertex *parameter, const double coefficient)
      : g2o::BaseUnaryEdge<2, Eigen::Vector2d, WayPointVertex>() {
    _vertices[0] = parameter;
    _information = coefficient * Eigen::Matrix2d::Identity();
  }
  void computeError() {
    _error = static_cast<const WayPointVertex *>(_vertices[0])->estimate() -
             _measurement;
  }
  void linearizeOplus() override {
    std::get<0>(_jacobianOplus) = Eigen::Matrix2d::Identity();
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class WayPointConstraint
    : public g2o::BaseUnaryEdge<1, double, WayPointVertex> {
public:
  const std::shared_ptr<std::vector<std::vector<Eigen::Vector2d>>> boundaries;
  const double agent_radius;
  WayPointConstraint(
      WayPointVertex *point, double _agent_radius,
      const std::shared_ptr<std::vector<std::vector<Eigen::Vector2d>>>
          &_boundaries,
      const double coefficient)
      : g2o::BaseUnaryEdge<1, double, WayPointVertex>(),
        boundaries(_boundaries), agent_radius(_agent_radius) {
    _vertices[0] = point;
    _information[0] = coefficient;
  }
  void setCoefficient(const double coefficient) {
    _information[0] = coefficient;
  }
  void computeError() {
    Eigen::Vector2d point =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    _error[0] = 0.;
    for (auto &boundary : *boundaries) {
      double distance = 2. * agent_radius;
      for (int i = 0; i < boundary.size(); i++) {
        Eigen::Vector2d point_a = boundary[i],
                        point_b = boundary[(i + 1) % boundary.size()];
        distance = std::min(distance, (point - point_a).norm());
        double length = (point_b - point_a).norm(),
               ip = (point_b - point_a).dot(point - point_a) / length;
        if (0 < ip && ip < length) {
          Eigen::Vector2d foot = point_a + (point_b - point_a) / length * ip;
          distance = std::min(distance, (point - foot).norm());
        }
      }
      _error[0] += distance < agent_radius
                       ? 1. / (distance + EPS) - 1. / agent_radius
                       : 0.;
    }
  }
  void linearizeOplus() override {
    Eigen::Vector2d point =
        static_cast<const WayPointVertex *>(_vertices[0])->estimate();
    std::get<0>(_jacobianOplus) = Eigen::MatrixXd::Zero(1, 2);
    for (auto &boundary : *boundaries) {
      double distance = 2. * agent_radius;
      Eigen::Vector2d nearest;
      for (int i = 0; i < boundary.size(); i++) {
        Eigen::Vector2d point_a = boundary[i],
                        point_b = boundary[(i + 1) % boundary.size()];
        if (distance > (point - point_a).norm()) {
          distance = (point - point_a).norm();
          nearest = point_a;
        }
        double length = (point_b - point_a).norm(),
               ip = (point_b - point_a).dot(point - point_a) / length;
        if (0 < ip && ip < length) {
          Eigen::Vector2d foot = point_a + (point_b - point_a) / length * ip;
          if (distance > (point - foot).norm()) {
            distance = (point - foot).norm();
            nearest = foot;
          }
        }
      }
      if (distance < agent_radius) {
        double dedd = -std::pow(distance + EPS, -2);
        std::get<0>(_jacobianOplus) +=
            (dedd * (point - nearest) / (distance + EPS)).transpose();
      }
    }
  }
  virtual bool read(std::istream &is) { return true; }
  virtual bool write(std::ostream &os) const { return true; }
};

class OurPreIteration : public g2o::HyperGraphAction {
public:
  const int pre_iteration, growth_iteration;
  const double growth_rate;
  const std::vector<CollisionConstraint *> c_constraints;

  OurPreIteration(const int _pre_iteration, const int _growth_iteration,
                  const double _growth_rate,
                  const std::vector<CollisionConstraint *> &_c_constraints)
      : g2o::HyperGraphAction(), pre_iteration(_pre_iteration),
        growth_iteration(_growth_iteration), growth_rate(_growth_rate),
        c_constraints(_c_constraints) {}

  g2o::HyperGraphAction *
  operator()(const g2o::HyperGraph *graph,
             g2o::HyperGraphAction::Parameters *parameters = 0) {
    assert(dynamic_cast<const g2o::SparseOptimizer *>(graph) &&
           "graph is not a SparseOptimizer");
    assert(dynamic_cast<g2o::HyperGraphAction::ParametersIteration *>(
               parameters) &&
           "error casting parameters");

    const g2o::SparseOptimizer *optimizer =
        static_cast<const g2o::SparseOptimizer *>(graph);
    g2o::HyperGraphAction::ParametersIteration *params =
        static_cast<g2o::HyperGraphAction::ParametersIteration *>(parameters);

    if (params->iteration >= pre_iteration &&
        params->iteration < pre_iteration + growth_iteration) {
      for (auto &constraint : c_constraints) {
        constraint->setInformation(growth_rate * constraint->information());
      }
    }

    return this;
  }
};

class OurPostIteration : public g2o::SparseOptimizerTerminateAction {
public:
  const int _minIteration;

  OurPostIteration(const int _min_iteration)
      : g2o::SparseOptimizerTerminateAction(), _minIteration(_min_iteration) {}

  g2o::HyperGraphAction *
  operator()(const g2o::HyperGraph *graph,
             g2o::HyperGraphAction::Parameters *parameters = 0) {
    using namespace g2o;

    assert(dynamic_cast<const SparseOptimizer *>(graph) &&
           "graph is not a SparseOptimizer");
    assert(dynamic_cast<HyperGraphAction::ParametersIteration *>(parameters) &&
           "error casting parameters");

    const SparseOptimizer *optimizer =
        static_cast<const SparseOptimizer *>(graph);
    HyperGraphAction::ParametersIteration *params =
        static_cast<HyperGraphAction::ParametersIteration *>(parameters);

    const_cast<SparseOptimizer *>(optimizer)->computeActiveErrors();
    if (params->iteration < 0) {
      // let the optimizer run for at least one iteration
      // Hence, we reset the stop flag
      setOptimizerStopFlag(optimizer, false);
    } else if (params->iteration < _minIteration) {
      // first iteration, just store the chi2 value
      _lastChi = optimizer->activeRobustChi2();
    } else {
      // compute the gain and stop the optimizer in case the
      // gain is below the threshold or we reached the max
      // number of iterations
      bool stopOptimizer = false;
      if (params->iteration < _maxIterations) {
        number_t currentChi = optimizer->activeRobustChi2();
        number_t gain = (_lastChi - currentChi) / currentChi;
        _lastChi = currentChi;
        if (gain >= 0 && gain < _gainThreshold)
          stopOptimizer = true;
      } else {
        stopOptimizer = true;
      }
      if (stopOptimizer) { // tell the optimizer to stop
        setOptimizerStopFlag(optimizer, true);
      }
    }
    return this;
  }
};

class TrajectoryOptimizer {
public:
  // Constants
  int number_of_way_points, number_of_agents,
      max_iteration = std::numeric_limits<int>::max();
  double agent_radius, max_velocity, max_acceleration;
  // bool consider_acceleration;
  // int path_begin, path_end;

  // Variables and constraints
  std::vector<std::vector<WayPointVertex *>> paths;
  ParameterVertex *delta;

  // Cost Parameters
  double collision_cost_const, velocity_cost_const, acceleration_cost_const,
      smooth_cost_const, smooth_acceleration_cost_const;

  g2o::SparseOptimizer optimizer;

  g2o::OptimizationAlgorithmLevenberg *solver;

  std::vector<PathCostConstraintMaxVelocity *> v_constraints;
  std::vector<PathCostConstraintMaxAcceleration *> a_constraints;
  std::vector<CollisionConstraint *> c_constraints;
  std::vector<WayPointConstraint *> w_constraints;

  int pre_iteration = 0, growth_iteration = 0;
  double pre_growth_rate = 1., growth_rate = 1., initial_lambda = 0.;

  bool variale_delta;

  double way_point_target_coefficient = 0.;
  std::vector<WayPointTargetConstraint *> start_constraints, goal_constraints;

  TrajectoryOptimizer(
      const YAML::Node &config, const int number_of_agents,
      const std::shared_ptr<std::vector<std::vector<Eigen::Vector2d>>>
          &boundaries) {
    number_of_way_points = config["number_of_way_points"].as<int>();
    this->number_of_agents = number_of_agents;
    if (config["max_iteration"]) {
      max_iteration = config["max_iteration"].as<int>();
    }

    typedef g2o::BlockSolver<g2o::BlockSolverTraits<-1, -1>> ABlockSolver;
    typedef g2o::LinearSolverEigen<ABlockSolver::PoseMatrixType> ALinearSolver;
    auto linearSolver = std::make_unique<ALinearSolver>();
    linearSolver->setBlockOrdering(false);
    solver = new g2o::OptimizationAlgorithmLevenberg(
        std::make_unique<ABlockSolver>(std::move(linearSolver)));
    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(config["verbose"] ? config["verbose"].as<bool>()
                                           : false);
    if (config["initial_lambda"]) {
      initial_lambda = config["initial_lambda"].as<double>();
      solver->setUserLambdaInit(initial_lambda);
    }

    v_constraints.clear();
    a_constraints.clear();
    c_constraints.clear();
    w_constraints.clear();

    int vertex_id = 0;
    delta = new ParameterVertex;
    delta->setId(vertex_id++);
    optimizer.addVertex(delta);
    variale_delta =
        config["variable_delta"] ? config["variable_delta"].as<bool>() : true;
    if (variale_delta) {
      auto delta_constraint = new ParameterTargetConstraint(
          delta, 0.0, std::pow((double)number_of_way_points, 2));
      optimizer.addEdge(delta_constraint);
    } else {
      delta->setFixed(true);
      delta->setEstimate(1. / number_of_way_points);
    }

    agent_radius = config["agent_radius"].as<double>();
    if (config["max_velocity"]) {
      max_velocity = config["max_velocity"].as<double>();
      velocity_cost_const = config["velocity_cost_const"].as<double>();
    } else {
      max_velocity = -1.;
    }
    if (config["max_acceleration"]) {
      max_acceleration = config["max_acceleration"].as<double>();
      acceleration_cost_const = config["acceleration_cost_const"].as<double>();
    } else {
      max_acceleration = -1.;
    }
    paths.resize(number_of_agents);
    // std::string path_cost_method =
    // config["path_cost_method"].as<std::string>();

    smooth_cost_const = config["smooth_cost_const"]
                            ? config["smooth_cost_const"].as<double>()
                            : -1.;
    smooth_acceleration_cost_const =
        config["smooth_acceleration_cost_const"]
            ? config["smooth_acceleration_cost_const"].as<double>()
            : -1.;

    /*
    consider_acceleration =
        (max_acceleration > 0. || smooth_acceleration_cost_const > 0.);
    if (!consider_acceleration) {
      path_begin = 0;
      path_end = number_of_way_points;
    } else {
      path_begin = 1;
      path_end = number_of_way_points + 1;
    }
    */

    collision_cost_const = config["collision_cost_const"].as<double>();

    for (int agent = 0; agent < number_of_agents; agent++) {
      auto &path = paths[agent];
      // int path_size = !consider_acceleration ? number_of_way_points + 1
      //                                       : number_of_way_points + 3;
      path.resize(number_of_way_points + 1);
      for (int i = 0; i <= number_of_way_points; i++) {
        path[i] = new WayPointVertex;
        path[i]->setId(vertex_id++);
        optimizer.addVertex(path[i]);
        auto w_constraint =
            new WayPointConstraint(path[i], agent_radius, boundaries,
                                   collision_cost_const / number_of_way_points);
        optimizer.addEdge(w_constraint);
        w_constraints.push_back(w_constraint);
      }
      // if (!consider_acceleration) {
      path[0]->setFixed(true);
      path[number_of_way_points]->setFixed(true);
      /*} else {
        path[0]->setFixed(true);
        path[1]->setFixed(true);
        path[number_of_way_points + 1]->setFixed(true);
        path[number_of_way_points + 2]->setFixed(true);
        }*/
      if (max_velocity > 0.) {
        for (int i = 0; i < number_of_way_points; i++) {
          auto path_cost_constraint = new PathCostConstraintMaxVelocity(
              path[i], path[i + 1], delta, max_velocity,
              velocity_cost_const / number_of_way_points);
          optimizer.addEdge(path_cost_constraint);
          v_constraints.push_back(path_cost_constraint);
        }
      }
      if (max_acceleration > 0.) {
        auto start_constraint = new PathCostConstraintMaxVelocity(
            path[0], path[1], delta, max_acceleration / number_of_way_points,
            acceleration_cost_const / (number_of_way_points + 1));
        optimizer.addEdge(start_constraint);
        v_constraints.push_back(start_constraint);
        for (int i = 1; i < number_of_way_points; i++) {
          auto path_cost_constraint = new PathCostConstraintMaxAcceleration(
              path[i - 1], path[i], path[i + 1], delta, max_acceleration,
              acceleration_cost_const / (number_of_way_points + 1));
          optimizer.addEdge(path_cost_constraint);
          a_constraints.push_back(path_cost_constraint);
        }
        auto goal_constraint = new PathCostConstraintMaxVelocity(
            path[number_of_way_points - 1], path[number_of_way_points], delta,
            max_acceleration / number_of_way_points,
            acceleration_cost_const / (number_of_way_points + 1));
        optimizer.addEdge(goal_constraint);
        v_constraints.push_back(goal_constraint);
      }

      if (smooth_cost_const > 0.) {
        for (int i = 0; i < number_of_way_points; i++) {
          auto path_cost_constraint = new PathCostConstraint(
              path[i], path[i + 1], delta,
              smooth_cost_const / 2 / number_of_way_points);
          optimizer.addEdge(path_cost_constraint);
        }
      }
      if (smooth_acceleration_cost_const > 0.) {
        auto start_constraint = new PathCostConstraint(
            path[0], path[1], delta,
            smooth_acceleration_cost_const / 2 * number_of_way_points);
        optimizer.addEdge(start_constraint);
        for (int i = 1; i < number_of_way_points; i++) {
          auto path_cost_constraint = new PathCostConstraintAcceleration(
              path[i - 1], path[i], path[i + 1], delta,
              smooth_acceleration_cost_const / 2 / number_of_way_points);
          optimizer.addEdge(path_cost_constraint);
        }
        auto goal_constraint = new PathCostConstraint(
            path[number_of_way_points - 1], path[number_of_way_points], delta,
            smooth_acceleration_cost_const / 2 * number_of_way_points);
        optimizer.addEdge(goal_constraint);
      }
    }
    for (int i = 0; i < number_of_agents; i++) {
      for (int j = i + 1; j < number_of_agents; j++) {
        for (int t = 1; t < number_of_way_points; t++) {
          auto collision_constraint = new CollisionConstraint(
              paths[i][t], paths[j][t], agent_radius,
              collision_cost_const / number_of_way_points);
          optimizer.addEdge(collision_constraint);
          c_constraints.push_back(collision_constraint);
        }
      }
    }
    if (config["pre_iteration"]) {
      pre_iteration = config["pre_iteration"].as<int>();
    }
    if (config["growth_iteration"]) {
      growth_iteration = config["growth_iteration"].as<int>();
      growth_rate = config["growth_rate"].as<double>();
    }
    auto pre_iteration_action = new OurPreIteration(
        pre_iteration, growth_iteration, growth_rate, c_constraints);
    optimizer.addPreIterationAction(pre_iteration_action);
    auto post_iteration_action =
        new OurPostIteration(pre_iteration + growth_iteration);
    if (config["terminal_threshold"]) {
      post_iteration_action->setGainThreshold(
          config["terminal_threshold"].as<double>());
    }
    optimizer.addPostIterationAction(post_iteration_action);

    if (config["way_point_target_coefficient"]) {
      way_point_target_coefficient =
          config["way_point_target_coefficient"].as<double>();
    }

    start_constraints.resize(number_of_agents);
    goal_constraints.resize(number_of_agents);
    for (int i = 0; i < number_of_agents; i++) {
      start_constraints[i] = new WayPointTargetConstraint(
          paths[i][0], way_point_target_coefficient);
      goal_constraints[i] = new WayPointTargetConstraint(
          paths[i][number_of_way_points], way_point_target_coefficient);
      optimizer.addEdge(start_constraints[i]);
      optimizer.addEdge(goal_constraints[i]);
    }
  }

  void growAccelerationCoefficients(const double rate) {
    for (auto &constraint : a_constraints) {
      constraint->setInformation(rate * constraint->information());
    }
  }

  void growCoefficients() {
    for (auto &constraint : v_constraints) {
      constraint->setInformation(growth_rate * constraint->information());
    }
    for (auto &constraint : a_constraints) {
      constraint->setInformation(growth_rate * constraint->information());
    }
    for (auto &constraint : c_constraints) {
      constraint->setInformation(growth_rate * constraint->information());
    }
    for (auto &constraint : w_constraints) {
      constraint->setInformation(growth_rate * constraint->information());
    }
  }

  void resetCoefficient() {
    for (auto &constraint : v_constraints) {
      constraint->setCoefficient(velocity_cost_const);
    }
    for (auto &constraint : a_constraints) {
      constraint->setCoefficient(acceleration_cost_const);
    }
    for (auto &constraint : c_constraints) {
      constraint->setCoefficient(collision_cost_const / number_of_way_points);
    }
    for (auto &constraint : w_constraints) {
      constraint->setCoefficient(collision_cost_const / number_of_way_points);
    }
  }

  void optimize(const std::vector<std::vector<Eigen::Vector2d>> &plan) {
    assert(plan.size() == number_of_agents);
    int makespan = plan[0].size() - 1;
    double initial_delta = (double)makespan / number_of_way_points;
    if (variale_delta) {
      delta->setEstimate(initial_delta);
    }
    for (int i = 0; i < number_of_agents; i++) {
      assert(plan[i].size() == makespan + 1);
      // if (!consider_acceleration) {
      paths[i][0]->setEstimate(plan[i][0]);
      paths[i][0]->setFixed(true);
      paths[i][number_of_way_points]->setEstimate(plan[i][makespan]);
      paths[i][number_of_way_points]->setFixed(true);
      /*} else {
      paths[i][0]->setEstimate(plan[i][0]);
      paths[i][1]->setEstimate(plan[i][0]);
      paths[i][number_of_way_points + 1]->setEstimate(plan[i][makespan]);
      paths[i][number_of_way_points + 2]->setEstimate(plan[i][makespan]);
      }*/
      for (int t = 1; t < number_of_way_points; t++) {
        double time = initial_delta * t;
        int s = std::floor(time);
        double r = time - s;
        paths[i][t]->setEstimate((1 - r) * plan[i][s] + r * plan[i][s + 1]);
      }
    }
    optimizer.initializeOptimization();
    resetCoefficient();
    int iterations = 0;
    while (iterations < max_iteration) {
      solver->setUserLambdaInit(initial_lambda);
      iterations += optimizer.optimize(max_iteration - iterations);
    }
  }

  std::vector<std::vector<Eigen::Vector2d>> getPlan() {
    std::vector<std::vector<Eigen::Vector2d>> plan(
        number_of_agents,
        std::vector<Eigen::Vector2d>(number_of_way_points + 1));
    for (int i = 0; i < number_of_agents; i++) {
      for (int t = 0; t <= number_of_way_points; t++) {
        plan[i][t] = paths[i][t]->estimate();
      }
    }
    return plan;
  }

  double getDelta() { return delta->estimate(); }

  void move_starts_and_goals(const std::vector<Eigen::Vector2d> &starts,
                             std::vector<Eigen::Vector2d> &goals) {
    for (int i = 0; i < number_of_agents; i++) {
      paths[i][0]->setFixed(false);
      start_constraints[i]->setMeasurement(starts[i]);
      paths[i][number_of_way_points]->setFixed(false);
      goal_constraints[i]->setMeasurement(goals[i]);
    }
    optimizer.initializeOptimization();
    resetCoefficient();
    optimizer.optimize(max_iteration);
  }

  double get_cost() {
    for (int i = 0; i < number_of_agents; i++) {
      paths[i][0]->setFixed(true);
      paths[i][number_of_way_points]->setFixed(true);
    }
    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();
    return optimizer.activeRobustChi2();
  }
};

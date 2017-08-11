#include <acado_gnuplot.hpp>
#include <acado_toolkit.hpp>
#include <iostream>

int main() {
  USING_NAMESPACE_ACADO;

  // VARIABLES:
  // -------------
  DifferentialState p;
  DifferentialState q;
  DifferentialState r;
  DifferentialState p_dot;
  DifferentialState q_dot;
  DifferentialState r_dot;

  Control db;
  Control dr1;
  Control dr2;
  Control dr3;
  Control dr4;

  Disturbance W;

  double Ts = 0.1;
  double u = 2.0;
  double delta_t = 0.1;
  double dup = -6.12;
  double dur = 0.375;
  double euq = -2.89;
  double fup = -0.245;
  double fur = -2.11;
  double ddr1 = 0.0478;
  double ddr2 = -0.0478;
  double ddr3 = 0.0478;
  double ddr4 = -0.0478;
  double edb = 0.0132;
  double edr1 = 0.0132;
  double edr2 = 0.0132;
  double edr3 = 0.0132;
  double edr4 = 0.0132;
  double fdr1 = -0.0132;
  double fdr2 = 0.0132;
  double fdr3 = 0.0132;
  double fdr4 = -0.0132;

  // DIFFERENTIAL EQUATION:
  // ----------------------
  DiscretizedDifferentialEquation f(Ts), fsim(Ts);

  // Controller inner model
  f << next(p) == p + p_dot * delta_t;
  f << next(q) == q + q_dot * delta_t;
  f << next(r) == r + r_dot * delta_t;
  f << next(p_dot) == dup * u / 100.0 * p + dur * u / 100.0 * r +
                          1.0 / 100.0 *
                              (ddr1 * u * u * dr1 + ddr2 * u * u * dr2 +
                               ddr3 * u * u * dr3 + ddr4 * u * u * dr4);
  f << next(q_dot) ==
      euq * u / 100.0 * q_dot +
          1.0 / 100.0 *
              (edb * u * u * db + edr1 * u * u * dr1 + edr2 * u * u * dr2 +
               edr3 * u * u * dr3 + edr4 * u * u * dr4);
  f << next(r_dot) == fup * u / 100.0 * p + fur * u / 100.0 * r +
                          1.0 / 100.0 *
                              (fdr1 * u * u * dr1 + fdr2 * u * u * dr2 +
                               fdr3 * u * u * dr3 + fdr4 * u * u * dr4);

  // Plant model
  fsim << next(p) == p + p_dot * delta_t + W * 1e-5;
  fsim << next(q) == q + q_dot * delta_t + W * 1e-5;
  fsim << next(r) == r + r_dot * delta_t + W * 1e-5;
  fsim << next(p_dot) == dup * u / 100.0 * p + dur * u / 100.0 * r +
                             1.0 / 100.0 *
                                 (ddr1 * u * u * dr1 + ddr2 * u * u * dr2 +
                                  ddr3 * u * u * dr3 + ddr4 * u * u * dr4);
  fsim << next(q_dot) ==
      euq * u / 100.0 * q_dot +
          1.0 / 100.0 *
              (edb * u * u * db + edr1 * u * u * dr1 + edr2 * u * u * dr2 +
               edr3 * u * u * dr3 + edr4 * u * u * dr4);
  fsim << next(r_dot) == fup * u / 100.0 * p + fur * u / 100.0 * r +
                             1.0 / 100.0 *
                                 (fdr1 * u * u * dr1 + fdr2 * u * u * dr2 +
                                  fdr3 * u * u * dr3 + fdr4 * u * u * dr4);

  // DEFINE LEAST SQUARE FUNCTION:
  // -----------------------------
  Function h;

  h << p;
  h << q;
  h << r;
  h << p_dot;
  h << q_dot;
  r << r_dot;
  h << db;
  h << dr1;
  h << dr2;
  h << dr3;
  h << dr4;

  DMatrix Q(11, 11);
  Q.setIdentity();

  DVector ref(11);

  // DEFINE AN OPTIMAL CONTROL PROBLEM:
  // ---------------------------------
  const double t_start = 0.0;
  const double t_end = 5.0;

  OCP ocp(t_start, t_end, 25);

  ocp.minimizeLSQ(Q, h, ref);
  ocp.subjectTo(f);
  //   ocp.subjectTo(-1 <= p <= 1);
  //   ocp.subjectTo(-1 <= q <= 1);
  //   ocp.subjectTo(-1 <= r <= 1);
  //   ocp.subjectTo(-1 <= p_dot <= 1);
  //   ocp.subjectTo(-1 <= q_dot <= 1);
  //   ocp.subjectTo(-1 <= r_dot <= 1);
  //   ocp.subjectTo(-1 <= db <= 1);
  //   ocp.subjectTo(-1 <= dr1 <= 1);
  //   ocp.subjectTo(-1 <= dr2 <= 1);
  //   ocp.subjectTo(-1 <= dr3 <= 1);
  //   ocp.subjectTo(-1 <= dr4 <= 1);

  // SETTIONG UP (SIMULATED) PROCESS:
  // --------------------------------
  OutputFcn identity;
  DynamicSystem dynamicSystem(fsim, identity);
  Process process(dynamicSystem, INT_RK45);

  VariablesGrid disturbance;
  disturbance.read("dist.txt");
  if (process.setProcessDisturbance(disturbance) != SUCCESSFUL_RETURN)
    exit(EXIT_FAILURE);

  // SETTING UP THE MPC CONTROLLER:
  // ------------------------------
  double samplingTime = 0.1;
  RealTimeAlgorithm alg(ocp, samplingTime);
  // alg.set(USE_REALTIME_ITERATIONS, NO);

  StaticReferenceTrajectory zeroReference;

  Controller controller(alg, zeroReference);

  DVector x0(6);
  x0.setZero();
  x0(0) = 1.0e-3;
  x0(1) = 1.0e-3;
  x0(2) = 1.0e-3;

  double startTime = 0.0;
  double endTime = 20.0;

  DVector uCon(5);
  uCon.setZero();
  VariablesGrid ySim;

  //   uint aaa = 0;
  //   aaa = controller.getNU();
  //   std::cout << aaa << std::endl;

  if (controller.init(startTime, x0) != SUCCESSFUL_RETURN)
    exit(EXIT_FAILURE);
  controller.getU(uCon);

  if (process.init(startTime, x0, uCon) != SUCCESSFUL_RETURN)
    exit(EXIT_FAILURE);
  process.getY(ySim);

  //   int nSteps = 0;
  //   double currentTime = startTime;

  //   while (currentTime <= endTime) {
  //     printf("\n*** Simulation Loop No. %d (starting at time %.3f) ***\n",
  //     nSteps,
  //            currentTime);
  //     if (controller.step(currentTime, ySim.getLastVector()) !=
  //     SUCCESSFUL_RETURN)
  //       exit(EXIT_FAILURE);
  //     controller.getU(uCon);
  //     std::cout << uCon;
  //     if (process.step(currentTime, currentTime + samplingTime, uCon) !=
  //         SUCCESSFUL_RETURN)
  //       exit(EXIT_FAILURE);
  //     process.getY(ySim);
  //     std::cout << ySim;

  //     ++nSteps;
  //     currentTime = (double)nSteps * samplingTime;
  //   }
  return EXIT_SUCCESS;
}
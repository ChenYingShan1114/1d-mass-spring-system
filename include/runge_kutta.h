//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    Eigen::VectorXd f1, f2, f3, f4;
    Eigen::VectorXd qdot1, qdot2, qdot3, qdot4;
    Eigen::VectorXd q1, q2, q3, q4;

    qdot1 = qdot;
    q1 = q;
    force(f1, q1, qdot1);

    qdot2 = qdot + 0.5 * dt * f1 / mass;
    q2 = q + 0.5 * dt * qdot1;
    force(f2, q2, qdot2);

    qdot3 = qdot + 0.5 * dt * f2 / mass;
    q3 = q + 0.5 * dt * qdot2;
    force(f3, q3, qdot3);

    qdot4 = qdot + dt * f3 / mass;
    q4 = q + dt * qdot3;
    force(f4, q4, qdot4);

    qdot = qdot + dt / 6 * (f1 + 2 * f2 + 2 * f3 + f4) / mass;
    q = q + dt / 6 * (qdot1 + 2 * qdot2 + 2 * qdot3 + qdot4);
}
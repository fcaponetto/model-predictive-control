#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "MPC.h"

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

int main() {
    uWS::Hub h;

    // MPC is initialized here!
    MPC mpc;

    h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                       uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);
        std::cout << sdata << std::endl;
        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {

                    ///**************************************************************
                    ///* GET THE CURRENT STATE
                    ///**************************************************************
                    // j[1] is the data JSON object
                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];
                    double px = j[1]["x"];
                    double py = j[1]["y"];
                    double psi = j[1]["psi"];
                    double v = j[1]["speed"];

                    ///**************************************************************
                    ///* CONVERT WAYPOINTS from GLOBAL SPACE to VEHICLE SPACE as VectorXd
                    ///**************************************************************
                    const int NUMBER_OF_WAYPOINTS = ptsx.size();
                    Eigen::VectorXd waypoints_xs(NUMBER_OF_WAYPOINTS);
                    Eigen::VectorXd waypoints_ys(NUMBER_OF_WAYPOINTS);

                    for (int i = 0; i < NUMBER_OF_WAYPOINTS; i++)
                    {
                        double dtx = ptsx[i] - px;
                        double dty = ptsy[i] - py;

                        waypoints_xs[i] = dtx * cos(psi) + dty * sin(psi);
                        waypoints_ys[i] = dty * cos(psi) - dtx * sin(psi);
                    }

                    ////*************************************************************
                    ///* FIT POLYNOMAL
                    ///**************************************************************
                    const int ORDER = 3;
                    auto coeffs = polyfit(waypoints_xs, waypoints_ys, ORDER);

                    ///**************************************************************
                    ///* GET POINTS TO DISPLAY FROM OUR FITTED POLYNOMIAL (ROAD CURVE)
                    ///**************************************************************
//                    std::vector<double> next_xs(N);
//                    std::vector<double> next_ys(N);
//                    const double D = 5.0;
//
//                    for (int i = 0; i < N; ++i) {
//
//                        const double dx = D * i;
//                        const double dy = K[3] * dx * dx * dx + K[2] * dx * dx + K[1] * dx + K[0];
//
//                        next_xs[i] = dx;
//                        next_ys[i] = dy;
//                    }

                    ///**************************************************************
                    ///* GENERATE CURRENT ERROR ESTIMATES (cte, epsi)
                    ///**************************************************************

                    // Estimate cross-track error (horizontal works reasonably well unless there's a lot of warpage
                    double cte = polyeval(coeffs, 0);
                    // Calculate orientation error
                    // TODO: check derivation is correct
                    double epsi = -atan(coeffs[1]);

                    ///**************************************************************
                    ///* GET THE CURRENT DELAYED STATE
                    ///* PREDICTION of the 100ms STATE into the future before it is sent to the solver.
                    ///**************************************************************

                    // Previous steering angle and throttle
                    double delta = j[1]["steering_angle"];
                    const double prev_a = j[1]["throttle"];

                    const double dt = 0.1;
                    const double Lf = 2.67;

                    // current state must be in vehicle coordinates with the delay factored in
                    // kinematic model is at play here
                    // note that at current state at vehicle coordinates:
                    // px, py, psi = 0.0, 0.0, 0.0
                    // note that in vehicle coordinates it is going straight ahead the x-axis
                    // which means position in vehicle's y-axis does not change
                    // the steering angle is negative the given value as we have
                    // as recall that during transformation we rotated all waypoints by -psi
                    const double current_px = 0.0 + v * dt;
                    const double current_py = 0.0;
                    const double current_psi = 0.0 + v * (-delta) / Lf * dt;
                    const double current_v = v + prev_a * dt;
                    const double current_cte = cte + v * sin(epsi) * dt;
                    const double current_epsi = epsi + v * (-delta) / Lf * dt;

                    const int NUMBER_OF_STATES = 6;
                    Eigen::VectorXd state(NUMBER_OF_STATES);
                    state <<
                        current_px,
                        current_py,
                        current_psi,
                        current_v,
                        current_cte,
                        current_epsi;

                    ///**************************************************************
                    ///* DETERMINE NEXT COURSE OF ACTION AND PREDICTED STATES
                    ///* USING MODEL PREDICTIVE CONTROL
                    ///**************************************************************

                    mpc.Solve(state, coeffs);

                    json msgJson;
                    // NOTE: Remember to divide by deg2rad(25) before you send the
                    //   steering value back. Otherwise the values will be in between
                    //   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
                    msgJson["steering_angle"] = mpc.steer;
                    msgJson["throttle"] = mpc.throttle;

                    // Display the MPC predicted trajectory
                    vector<double> mpc_x_vals;
                    vector<double> mpc_y_vals;

                    /**
                     * TODO: add (x,y) points to list here, points are in reference to
                     *   the vehicle's coordinate system the points in the simulator are
                     *   connected by a Green line
                     */

                    msgJson["mpc_x"] = mpc.future_xs;
                    msgJson["mpc_y"] = mpc.future_ys;

                    // Display the waypoints/reference line
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    /**
                     * TODO: add (x,y) points to list here, points are in reference to
                     *   the vehicle's coordinate system the points in the simulator are
                     *   connected by a Yellow line
                     */

                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;


                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                    std::cout << msg << std::endl;
                    // Latency
                    // The purpose is to mimic real driving conditions where
                    //   the car does actuate the commands instantly.
                    //
                    // Feel free to play around with this value but should be to drive
                    //   around the track with 100ms latency.
                    //
                    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }  // end "telemetry" if
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }  // end websocket if
    }); // end h.onMessage

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }

    h.run();
}
/****************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC                                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of 3dThesis. 3dThesis is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <cmath>
#include <algorithm>
#include "DataStructs.h"
#include "Calc.h"
#include "Util.h"

#include <iostream>
#include <fstream>

void Calc::Integrate_Parallel(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	const int numThreads = sim.settings.thnum;
	static vector<Nodes> nodes_par;
	if (numThreads>1) {
		if (!nodes_par.size()) {
			nodes_par.resize(numThreads);
			#pragma omp parallel num_threads(numThreads)
			{
				Nodes th_nodes;
				#pragma omp for schedule(static)
				for (int i = 0; i < numThreads; i++) {
					Calc::Integrate_Serial(th_nodes, sim, t+i*sim.param.dt, isSol);
					nodes_par[numThreads - i - 1] = th_nodes;
				}
			}
		}
		nodes = nodes_par.back();
		nodes_par.pop_back();
	}
	else {
		Calc::Integrate_Serial(nodes, sim, t, isSol);
	}
	return;
}

void Calc::Integrate_Serial(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	if (sim.settings.compress) { 
		Calc::GaussCompressIntegrate(nodes, sim, t, isSol); 
		if (!isSol) {
			Nodes nodes_reg;
			Calc::GaussIntegrate(nodes_reg, sim, t, isSol);
			if (nodes_reg.size <= nodes.size) { nodes = nodes_reg; }
		}
	}
	else { 
		Calc::GaussIntegrate(nodes, sim, t, isSol);
	}

	if (sim.domain.use_BCs) { Calc::AddBCs(nodes, sim.domain); }

	return;
}

void Calc::GaussIntegrate(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	
	// Cache the starting segment for the next integration call.
	static vector<int> start_seg(sim.paths.size(), 1);

	// Quadrature node locations for order 2, 4, 8, and 16
	static const double locs[30] = {
		-0.57735027,  0.57735027,
		-0.86113631, -0.33998104,  0.33998104,  0.86113631,
		-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
		-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
	};

	// Quadrature weights for order 2, 4, 8, and 16
	static const double weights[30] = {
		1.0, 1.0,
		0.34785485, 0.65214515, 0.65214515, 0.34785485,
		0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
		0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
	};
	for (int i = 0; i < sim.paths.size(); i++) {

		const Beam& beam = sim.beams[i];
		const vector<path_seg>& path = sim.paths[i];
		int seg_temp = start_seg[i];

		const double beta = pow(3.0 / PI, 1.5) * beam.q / (sim.material.rho * sim.material.cps);

		while ((t > path[seg_temp].seg_time) && (seg_temp + 1 < path.size())) { seg_temp++; }
		while ((t < path[seg_temp - 1].seg_time) && (seg_temp - 1 > 0)) { seg_temp--; }
		
		if (!isSol) {start_seg[i] = seg_temp;}

		const double t0 = Util::t0calc(t, beam, sim.material, sim.settings);
		
		double curStep_max_start = beam.nond_dt;
		
		if (isSol) { curStep_max_start *= beam.az; }

		double curStep_max = curStep_max_start;

		double curStep_use = curStep_max;
		
		int curOrder = 16;

		// Preserve the instantaneous heat-source contribution to the Laplacian.
		int_seg current_beam = Util::GetBeamLoc(t, seg_temp, path, sim); 
		current_beam.phix = (beam.ax * beam.ax + 0.0);
		current_beam.phiy = (beam.ay * beam.ay + 0.0);
		current_beam.phiz = (beam.az * beam.az + 0.0);
		current_beam.qmod *= beta;
		current_beam.dtau = 0.0;
		if (t <= path.back().seg_time && current_beam.qmod > 0.0) { Util::AddToNodes(nodes, current_beam); }

		bool tflag = true;
		double t2 = t;
		double t1 = t;
		double tpp = 0.0;

		if (t > path.back().seg_time) {
			tpp += t - path.back().seg_time;
			t2 = path.back().seg_time;
		}

		while (tpp >= 2 * curStep_max - curStep_max_start) {
			curStep_max *= 2.0;
			if (curOrder != 2) {
				curOrder = (curOrder / 2);
			}
		}

		while (tflag) {
			bool switchSeg = false;

			const double ref_time = Util::GetRefTime(tpp, seg_temp, path, beam);

			curStep_use = (ref_time < curStep_max) ? ref_time : curStep_max;

			t1 = t2 - curStep_use;

			const double next_time = path[seg_temp - 1].seg_time;

			if (t1 <= next_time) {
				t1 = next_time;
				if (next_time > t0) { switchSeg = true; }
				else { tflag = false; }
			}

			double tau, ct;
			for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
				double tp = 0.5 * ((t2 - t1) * locs[a] + (t2 + t1));
				tau = t - tp;
				ct = 12.0 * sim.material.a * tau;

				current_beam = Util::GetBeamLoc(tp, seg_temp, path, sim);
				current_beam.phix = (beam.ax * beam.ax + ct);
				current_beam.phiy = (beam.ay * beam.ay + ct);
				current_beam.phiz = (beam.az * beam.az + ct);
				current_beam.qmod *= beta;
				current_beam.dtau = 0.5 * (t2 - t1) * weights[a];

				if (current_beam.qmod > 0.0 && current_beam.dtau > 0.0) { Util::AddToNodes(nodes, current_beam); }
			}

			if (switchSeg) {seg_temp--;}
			
			tpp += (t2 - t1);

			t2 = t1;

			if (tpp >= 2 * curStep_max - curStep_max_start) {
				curStep_max *= 2.0;
				if (curOrder != 2) {
					curOrder = (curOrder / 2);
				}
			}
		}
	}
	return;
}

void Calc::GaussCompressIntegrate(Nodes& nodes, const Simdat& sim, const double t, const bool isSol) {
	// Cache the starting segment for the next integration call.
	static vector<int> start_seg(sim.paths.size(), 1);

	// Quadrature node locations for order 2, 4, 8, and 16
	static const double locs[30] = {
		-0.57735027,  0.57735027,
		-0.86113631, -0.33998104,  0.33998104,  0.86113631,
		-0.96028986, -0.79666648, -0.52553241, -0.18343464,  0.18343464,  0.52553241, 0.79666648,  0.96028986,
		-0.98940093, -0.94457502, -0.8656312, -0.75540441, -0.61787624, -0.45801678, -0.28160355, -0.09501251, 0.09501251, 0.28160355, 0.45801678, 0.61787624, 0.75540441, 0.8656312, 0.94457502, 0.98940093
	};

	// Quadrature weights for order 2, 4, 8, and 16
	static const double weights[30] = {
		1.0, 1.0,
		0.34785485, 0.65214515, 0.65214515, 0.34785485,
		0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854,
		0.02715246, 0.06225352, 0.09515851, 0.12462897, 0.14959599, 0.16915652,0.18260342, 0.18945061, 0.18945061, 0.18260342, 0.16915652, 0.14959599,0.12462897, 0.09515851, 0.06225352, 0.02715246
	};

	for (int i = 0; i < sim.paths.size(); i++) {

		const Beam& beam = sim.beams[i];
		const vector<path_seg>& path = sim.paths[i];
		int seg_temp = start_seg[i];

		const double beta = pow(3.0 / PI, 1.5) * beam.q / (sim.material.rho * sim.material.cps);

		while ((t > path[seg_temp].seg_time) && (seg_temp + 1 < path.size())) { seg_temp++; }
		while ((t < path[seg_temp - 1].seg_time) && (seg_temp - 1 > 0)) { seg_temp--; }

		if (!isSol) { start_seg[i] = seg_temp; }

		const double t0 = Util::t0calc(t, beam, sim.material, sim.settings);

		double curStep_max_start = beam.nond_dt;

		if (isSol) { curStep_max_start *= beam.az; }

		double curStep_max = curStep_max_start;

		double curStep_use = curStep_max;

		int curOrder = 16;

		// Preserve the instantaneous heat-source contribution to the Laplacian.
		int_seg current_beam = Util::GetBeamLoc(t, seg_temp, path, sim);
		current_beam.phix = (beam.ax * beam.ax + 0.0);
		current_beam.phiy = (beam.ay * beam.ay + 0.0);
		current_beam.phiz = (beam.az * beam.az + 0.0);
		current_beam.qmod *= beta;
		current_beam.dtau = 0.0;
		if (t <= path.back().seg_time && current_beam.qmod>0.0) { Util::AddToNodes(nodes, current_beam); }

		bool tflag = true;
		double t2 = t;
		double t1 = t;
		double tpp = 0.0;

		if (t > path.back().seg_time) {
			tpp += t - path.back().seg_time;
			t2 = path.back().seg_time;
		}

		while (tpp >= 2 * curStep_max - curStep_max_start) {
			curStep_max *= 2.0;
			if (curOrder != 2) {
				curOrder = (curOrder / 2);
			}
		}

		while (tflag) {

			double r2, dist2, xp, yp, xs, ys, dx, dy, ts, dt;
			double sum_t = 0, sum_qmodt = 0, sum_qmodtx = 0, sum_qmodty = 0;
			int seg_temp_2;
			int num_comb_segs = 0;

			int quit = 0;
			while (true) {
				xs = path[seg_temp].sx;
				ys = path[seg_temp].sy;
				ts = path[seg_temp].seg_time;
				if (ts <= t0) { quit = 1; break; }
				if (Util::InRMax(xs, ys, sim.domain, sim.settings)) { break; }
				else { tpp += t2 - ts; t2 = ts; }
			}
			if (quit) { break; }

			while (tpp >= 2 * curStep_max - curStep_max_start) {
				curStep_max *= 2.0;
				if (curOrder != 2) {
					curOrder = (curOrder / 2);
				}
			}

			double ref_time = Util::GetRefTime(tpp, seg_temp, path, beam);
			if (ref_time < curStep_max) { curStep_use = ref_time; }
			else { curStep_use = curStep_max; }

			int_seg current_beam_t2 = Util::GetBeamLoc(t2, seg_temp, path, sim);
			xp = current_beam_t2.xb;
			yp = current_beam_t2.yb;

			r2 = log(2.0) / 8.0 * (beam.ax * beam.ax) * (12.0 * (t - t2) * sim.material.a / (beam.ax * beam.ax) + 1.0);

			seg_temp_2 = seg_temp;

			t1 = path[seg_temp - 1].seg_time;

			if (t1 < t2 - curStep_use) {
				num_comb_segs = 0;
				t1 = t2 - curStep_use;
			}
			else {
				bool cflag = true;
				while (cflag) {
					if (path[seg_temp_2 - 1].seg_time <= t0) { tflag = 0; break; }

					xs = path[seg_temp_2 - 1].sx;
					ys = path[seg_temp_2 - 1].sy;
					dist2 = (xs - xp) * (xs - xp) + (ys - yp) * (ys - yp);
					ts = path[seg_temp_2 - 1].seg_time;

					if (!Util::InRMax(xp, yp, sim.domain, sim.settings)) { break; }

					// Stop compression beyond the diffusion distance or across a
					// power discontinuity; otherwise include the neighboring segment.
					if ((dist2 > r2) || (path[seg_temp_2].sqmod != path[seg_temp_2 - 1].sqmod && (curOrder > 2 || (t2 - ts) > (curStep_max / 64.0)))) {
						if (path[seg_temp_2].smode) {
							sum_qmodtx += path[seg_temp_2].sx * path[seg_temp_2].sqmod * (t1 - ts);
							sum_qmodty += path[seg_temp_2].sy * path[seg_temp_2].sqmod * (t1 - ts);
							sum_qmodt += path[seg_temp_2].sqmod * (t1 - ts);
							sum_t += (t1 - ts);
							t1 = ts;
						}
						else {
							dx = path[seg_temp_2].sx - path[seg_temp_2 - 1].sx;
							dy = path[seg_temp_2].sy - path[seg_temp_2 - 1].sy;
							dt = path[seg_temp_2].seg_time - ts;

							double t_int = ts + dt * (sqrt(((xp - xs) * (xp - xs) + (yp - ys) * (yp - ys)) / (dx * dx + dy * dy)) - sqrt(r2 / (dx * dx + dy * dy)));
							if (t_int < t2 - curStep_max) { t_int = t2 - curStep_max; }
							sum_qmodtx += (xs + dx * ((t1 + t_int) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * (t1 - t_int);
							sum_qmodty += (ys + dy * ((t1 + t_int) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * (t1 - t_int);
							sum_qmodt += path[seg_temp_2].sqmod * (t1 - t_int);
							sum_t += (t1 - t_int);
							t1 = t_int;
						}
						cflag = false;
					}
					else {
						if (ts < t2 - curStep_max) {
							cflag = false;
							if (path[seg_temp_2].smode) { t1 = t2 - curStep_max; }
							else { t1 = t2 - curStep_max; }
						}
						else {
							t1 = ts;
							num_comb_segs++;
						}

						if (path[seg_temp_2].smode) {
							dt = (path[seg_temp_2].seg_time - t1);

							sum_qmodtx += path[seg_temp_2].sx * path[seg_temp_2].sqmod * dt;
							sum_qmodty += path[seg_temp_2].sy * path[seg_temp_2].sqmod * dt;
							sum_qmodt += path[seg_temp_2].sqmod * dt;
							sum_t += dt;
						}
						else {
							dx = path[seg_temp_2].sx - path[seg_temp_2 - 1].sx;
							dy = path[seg_temp_2].sy - path[seg_temp_2 - 1].sy;
							dt = (path[seg_temp_2].seg_time - t1);

							sum_qmodtx += (xs + dx * ((path[seg_temp_2].seg_time + t1) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * dt;
							sum_qmodty += (ys + dy * ((path[seg_temp_2].seg_time + t1) / 2.0 - ts) / dt) * path[seg_temp_2].sqmod * dt;
							sum_qmodt += path[seg_temp_2].sqmod * dt;
							sum_t += dt;
						}
					}
					if (t1 == path[seg_temp_2 - 1].seg_time) {
						xs = path[seg_temp_2 - 1].sx;
						ys = path[seg_temp_2 - 1].sy;
						seg_temp_2--;
						if (!Util::InRMax(xs, ys, sim.domain, sim.settings)) { break; }
					}
				}

				if (!num_comb_segs) {
					t1 = path[seg_temp - 1].seg_time;
					seg_temp_2 = seg_temp - 1;
				}
			}

			double tau, ct;
			if (num_comb_segs) {
				int_seg current_beam;
				if (sum_qmodt > 0.0) {
					current_beam.xb = sum_qmodtx / sum_qmodt;
					current_beam.yb = sum_qmodty / sum_qmodt;
					current_beam.qmod = sum_qmodt / sum_t;		
					for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
						double tp = 0.5 * ((t2 - t1) * locs[a] + (t2 + t1));
						tau = t - tp;
						ct = 12.0 * sim.material.a * tau;

						current_beam.phix = (beam.ax * beam.ax + ct);
						current_beam.phiy = (beam.ay * beam.ay + ct);
						current_beam.phiz = (beam.az * beam.az + ct);
						current_beam.qmod *= beta;
						current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
						if ((current_beam.qmod > 0.0) && (current_beam.dtau > 0.0)) { Util::AddToNodes(nodes, current_beam); }
					}
				}

			}
			else {
				for (int a = (2 * curOrder - 3); a > (curOrder - 3); a--) {
					double tp = 0.5 * ((t2 - t1) * locs[a] + (t2 + t1));	
					tau = t - tp;
					ct = 12.0 * sim.material.a * tau;

					int_seg current_beam = Util::GetBeamLoc(tp, seg_temp, path, sim);
					current_beam.phix = (beam.ax * beam.ax + ct);
					current_beam.phiy = (beam.ay * beam.ay + ct);
					current_beam.phiz = (beam.az * beam.az + ct);
					current_beam.qmod *= beta;
					current_beam.dtau = 0.5 * (t2 - t1) * weights[a];
					if ((current_beam.qmod > 0.0) && (current_beam.dtau > 0.0)) { Util::AddToNodes(nodes, current_beam); }
				}
			}

			tpp += t2 - t1;
			t2 = t1;

			if (tpp >= 2 * curStep_max - curStep_max_start) {
				curStep_max *= 2.0;
				if (curOrder != 2) {
					curOrder = (curOrder / 2);
				}
			}

			seg_temp = seg_temp_2;
			if (t1 <= t0) { tflag = 0; }
		}
	}
	return;
}

void Calc::AddBCs(Nodes& nodes, const Domain& domain) {

	vector<vector<int>> allCoords;
	vector<vector<int>> newCoords;
	vector<int> coords = { 0,0,0 };
	newCoords.push_back(coords);
	allCoords.push_back(coords);

	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;

	xmin = domain.BC_xmin;
	xmax = domain.BC_xmax;
	ymin = domain.BC_ymin;
	ymax = domain.BC_ymax;
	zmin = domain.BC_zmin;
	zmax = domain.zmax;

	vector<vector<int>> checkCoords = newCoords;
	newCoords.clear();
	// Strength controls the number of image-source reflections.
	int iter = 0; int strength = 0; 
	if (domain.BC_reflections != INT_MAX) { strength = domain.BC_reflections; }
	while (checkCoords.size() > 0 && iter < strength) {
		for (int i = 0; i < checkCoords.size(); i++) {
			if (xmin != DBL_MAX) {
				coords = checkCoords[i];
				coords[0] = (-1 - coords[0]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {;
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (xmax != DBL_MAX) {
				coords = checkCoords[i];
				coords[0] = (1 - coords[0]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (ymin != DBL_MAX) {
				coords = checkCoords[i];
				coords[1] = (-1 - coords[1]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (ymax != DBL_MAX) {
				coords = checkCoords[i];
				coords[1] = (1 - coords[1]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (zmin != DBL_MAX) {
				coords = checkCoords[i];
				coords[2] = (-1 - coords[2]);
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end()) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
			if (zmax != DBL_MAX) {
				coords = checkCoords[i];
				coords[2] = (1 - coords[2]);
				
				if (std::find(allCoords.begin(), allCoords.end(), coords) == allCoords.end() && (coords[2]!=1)) {
					allCoords.push_back(coords);
					newCoords.push_back(coords);
				}
			};
		}
		checkCoords = newCoords;
		newCoords.clear();
		iter += 1;
	}

	int org_size = nodes.size;
	Nodes nodes_org = nodes;
	Util::ClearNodes(nodes);

	// Convert the reference coordinates back to physical coordinates.
	for (int i = 0; i < allCoords.size(); i++) {
		int xRef = allCoords[i][0];
		int yRef = allCoords[i][1];
		int zRef = allCoords[i][2];	

		int nX = std::abs(xRef % 2);
		int nY = std::abs(yRef % 2);
		int nZ = std::abs(zRef % 2);	

		Nodes nodes2 = nodes_org;
		for (int i = 0; i < org_size; i++) { 
			nodes2.xb[i] = (xRef + nX) * xmax - (xRef - nX) * xmin + (1 - 2 * nX) * nodes2.xb[i];
			nodes2.yb[i] = (yRef + nY) * ymax - (yRef - nY) * ymin + (1 - 2 * nY) * nodes2.yb[i];
			nodes2.zb[i] = (zRef + nZ) * zmax - (zRef - nZ) * zmin + (1 - 2 * nZ) * nodes2.zb[i];
		}
		Util::CombineNodes(nodes, nodes2);	
	}

	return;
}

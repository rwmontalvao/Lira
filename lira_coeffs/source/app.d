// Copyright 2021 Astex Therapeutics Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Code libraries
import geotools;

// Standard libraries
import std.conv;
import std.math;
import std.array;
import std.stdio;
import std.getopt;

// Linear Algebra
import gl3n.linalg;

alias vec3d = Vector!(double, 3);

int main(string[] args)
{
	if (args.length < 2)
	{
		writeln("Usage: lira_coeffs --in_file <file.tmesh> --out_file <file.coeffs> --rmsd");
		return 1;
	}
	// Maximum l value
	const int lmax = 16;

	// Number of points for
	// the Spherical Sampling
	const int n_samples = 10_000;

	//  Compute RMSD?
	bool compute_rmsd = false;

	//  Verbose?
	bool verbose = false;

	// Mesh in the tmesh format
	auto in_file_name = "";
	auto out_file_name = "";

	auto cmds = getopt(
		args,
		"in_file", &in_file_name,
		"out_file", &out_file_name,
		"rmsd", &compute_rmsd,
		"verbose", &verbose);

	if (cmds.helpWanted)
	{
		defaultGetoptPrinter("Lira 1.0 RID - Copyright (c) 2021 Rinaldo Wander Montalvão, PhD", cmds
				.options);
		return 0;
	}

	if (verbose)
	{
		writeln("Lira 1.0 RID - Copyright (c) 2021 Rinaldo Wander Montalvão, PhD");
		writeln("   Rega Institute for Medical Research - KU Leuven - Belgium");
		writeln("Synthetic Biology Through Directed Evolution - The Pinheiro Lab\n");
	}

	vec3d[] vertices;
	double[] charges;

	// Parse mesh vertices
	string[] lines;
	try
	{
		foreach (line; File(in_file_name).byLine())
		{
			immutable string record = to!string(line);

			lines ~= record;

			immutable auto tokens = record.split();
			if (tokens.length == 9)
			{
				vertices ~= vec3d(to!(double[3])(tokens[0 .. 3]));
				charges ~= to!double(tokens[7]);
			}
		}
	}
	catch (Exception exc)
	{
		stderr.writefln("Error: %s", exc.msg);
		return 1;
	}

	// Parse the triangle format
	int[3][] triangles;
	foreach (i; 0 .. lines.length)
	{
		immutable auto line = lines[i];
		// If the line marks a triangle
		// get the next 3 lines
		if (line == "3 0")
		{
			triangles ~= to!(int[3])(lines[i + 1 .. i + 4]);
		}
	}

	// Generate the spherical sampling (radius = 1)
	auto samples_rtp = spherical_samples(n_samples);

	/*
		Surface sampling and Ylm pre-computing
	*/
	// gFile out_file = File("surface.txt", "w");

	vec3d[] surface_rtp;
	double[] surface_charge;
	double[][] ylm_values;

	// Find the triangle that intersects each
	// spherical sample
	vec3d origin = vec3d(0.0, 0.0, 0.0);
	vec3d intersection = vec3d(0.0, 0.0, 0.0);
	foreach (vec3d coord_rtp; samples_rtp)
	{
		vec3d coord = rtp2xyz(coord_rtp);

		// Check each triangle
		foreach (int[3] triangle; triangles)
		{
			// Find triangle vertices elements and
			// their charges
			vec3d[3] elements;
			double[3] values;

			elements[0] = vertices[triangle[0]];
			elements[1] = vertices[triangle[1]];
			elements[2] = vertices[triangle[2]];

			values[0] = charges[triangle[0]];
			values[1] = charges[triangle[1]];
			values[2] = charges[triangle[2]];

			// Find ray/triangle intersecion
			if (ray_intersects_triangle(origin, coord, elements, intersection))
			{
				// Interpolate the charge value at the intersection
				immutable double[2] weights = barycentric_coordinates(elements, intersection);

				immutable double u = weights[0];
				immutable double v = weights[1];

				immutable double value = u * values[1] + v * values[2] + (1 - u - v) * values[0];

				// immutable vec3d test = u * elements[1] + v * elements[2] + (1 - u - v) * elements[0];
				// writeln(elements);
				// writeln(values);
				// writeln(intersection);
				// writefln("%f %f %f", u, v, value);
				// writeln(test);

				// out_file.writefln("%f, %f, %f, %f", intersection.x,	intersection.y, intersection.z, value);

				vec3d intersection_rtp = xyz2rtp(intersection);

				surface_rtp ~= intersection_rtp;
				surface_charge ~= value;

				// Compute the Spherical Harmonics for
				// the intersection vector
				double[] band;
				foreach (l; 0 .. lmax + 1)
				{
					foreach (m; -l .. l + 1)
					{
						band ~= ylm(l, m, intersection_rtp.y, intersection_rtp.z);
					}
				}
				ylm_values ~= band;

				break; // only one triangle intersects
			}
		}
	}
	// out_file.close();

	/**
	 * Compute the Spherical Harmonics Coefficients
	 */
	immutable double delta = 4.0 * PI / to!double(n_samples);

	double[] coeffs_surface, coeffs_surface_charge;

	int i = 0;
	foreach (l; 0 .. lmax + 1)
	{
		foreach (m; -l .. l + 1)
		{
			double sum_surface = 0.0;
			double sum_surface_charge = 0.0;
			foreach (j; 0 .. surface_rtp.length)
			{
				sum_surface += surface_rtp[j].x * ylm_values[j][i] * delta;
				sum_surface_charge += surface_charge[j] * ylm_values[j][i] * delta;
			}
			coeffs_surface ~= sum_surface;
			coeffs_surface_charge ~= sum_surface_charge;
			i++;
		}
	}

	/**
	 *	Compute the RIDs
	 */
	double[] rid_surface, rid_charge;
	i = 0;
	foreach (l; 0 .. lmax + 1)
	{
		double sum_surface = 0.0, sum_charge = 0.0;
		foreach (m; -l .. l + 1)
		{
			sum_surface += pow(coeffs_surface[i], 2);
			sum_charge += pow(coeffs_surface_charge[i], 2);
			i++;
		}
		rid_surface ~= sum_surface;
		rid_charge ~= sum_charge;
	}

	/** 
	 * Output the coeffs
	 */

	File rid_file = File(out_file_name, "w");

	rid_file.writeln(in_file_name);
	rid_file.write("sr:");
	rid_file.writeln(rid_surface);
	rid_file.write("cr:");
	rid_file.writeln(rid_charge);
	rid_file.write("sc:");
	rid_file.writeln(coeffs_surface);
	rid_file.write("cc:");
	rid_file.writeln(coeffs_surface_charge);
	rid_file.close();

	/**
	 *	Compute the RMSD for the reconstruction
	 */

	if (compute_rmsd)
	{
		File csv_file = File("rmsd_data.csv", "w");

		csv_file.writeln("radius,radius_sh,charge,charge_sh");

		double max_surface = surface_rtp[0].x;
		double min_surface = surface_rtp[0].x;

		double max_charge = surface_charge[0];
		double min_charge = surface_charge[0];

		double total_surface = 0.0, total_charge = 0.0;
		foreach (j; 0 .. surface_rtp.length)
		{

			immutable double radius = surface_rtp[j].x;
			// immutable double theta = surface_rtp[j].y;
			// immutable double phi = surface_rtp[j].z;

			if (radius > max_surface)
				max_surface = radius;
			if (radius < min_surface)
				min_surface = radius;

			immutable double charge = surface_charge[j];

			if (charge > max_charge)
				max_charge = charge;
			if (charge < min_charge)
				min_charge = charge;

			i = 0;
			double sum_surface = 0.0, sum_charge = 0.0;
			foreach (l; 0 .. lmax + 1)
			{
				foreach (m; -l .. l + 1)
				{
					sum_surface += coeffs_surface[i] * ylm_values[j][i];

					sum_charge += coeffs_surface_charge[i] * ylm_values[j][i];
					i++;
				}
			}
			//writefln("%f - %f", radius, sum_surface);
			total_surface += pow(radius - sum_surface, 2);
			total_charge += pow(charge - sum_charge, 2);
			csv_file.writefln("%f,%f,%f,%f", radius, sum_surface, charge, sum_charge);
		}

		csv_file.close();

		immutable double rmsd_surface = sqrt(total_surface / to!double(surface_rtp.length));
		immutable double rmsd_charge = sqrt(total_charge / to!double(surface_rtp.length));

		writefln("Min. = %f", min_surface);
		writefln("Max. = %f", max_surface);
		writefln("RMSD = %f\n", rmsd_surface);

		writefln("Min. = %f", min_charge);
		writefln("Max. = %f", max_charge);
		writefln("RMSD = %f", rmsd_charge);

	}

	return 0;
}

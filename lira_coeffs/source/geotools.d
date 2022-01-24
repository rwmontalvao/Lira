// Copyright (c) 2021, KU Leuven
// Licensed under the Non-Profit Open Software License version 3.0.
// SPDX-License-Identifier: NPOSL-3.0

module geotools;

import std.conv;
import std.math;
import std.random;

import gl3n.linalg;

alias vec3d = Vector!(double, 3);


/******************************************************
 * Spherical Samples for the unit Sphere
 * Params:
 *      n_samples = number of points
 * Returns:
 *      array of spherical coordinates (r, theta, phi)
 *******************************************************/
pure vec3d[] spherical_samples(in int n_samples)
{
    vec3d[] coords;

    immutable int n = to!int(sqrt(to!double(n_samples)));

    immutable double oneovern = 1.0 / to!double(n);

    // auto rnd = Mt19937(unpredictableSeed);
    auto rnd = Mt19937(42);

    foreach (a; 0 .. n)
    {
        foreach (b; 0 .. n)
        {
            immutable double x = (to!double(a) + uniform01(rnd)) * oneovern;
            immutable double y = (to!double(b) + uniform01(rnd)) * oneovern;

            immutable double theta = 2.0 * acos(sqrt(1.0 - x));
            immutable double phi = 2.0 * PI * y;

            coords ~= vec3d(1.0, theta, phi);
        }
    }

    return coords;
}


/************************************************************
 * Convert Spherical Coordinates to Rectangular Coordinates
 * Params:
 *      vec = array of spherical coordinates (r, theta, phi)
 * Returns:
 *      array of rectangular coordinates (z, y, z)
 ************************************************************/
pure nothrow vec3d rtp2xyz(in vec3d vec)
{
    immutable double x = vec.x * sin(vec.y) * cos(vec.z);
    immutable double y = vec.x * sin(vec.y) * sin(vec.z);
    immutable double z = vec.x * cos(vec.y);

    return vec3d(x, y, z);
}

/************************************************************
 * Convert Rectangular Coordinates to Spherical Coordinates
 * Params:
 *      vec = array of rectangular coordinates (z, y, z)
 * Returns:
 *      array of spherical coordinates (r, theta, phi)
 *************************************************************/
pure nothrow vec3d xyz2rtp(in vec3d vec)
{

    immutable double r = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    immutable double theta = acos(vec.z / r);
    immutable double phi = atan2(vec.y, vec.x);

    return vec3d(r, theta, phi);
}


/************************************************************
 * Test the intersection between a ray and a triangle
 * Params:
 *      rayOrigin = origin 3D vector
 *      rayVector = ray 3D vector 
 *      inTriangle = triangles 3D vectors
 *      outIntersectionPoint = intersection point 3D vector (changes)
 * Returns:
 *       bool test result
 *************************************************************/
pure nothrow bool ray_intersects_triangle(in vec3d rayOrigin, in vec3d rayVector,
        in vec3d[3] inTriangle, ref vec3d outIntersectionPoint)
{
    const double EPSILON = 0.0000001;

    immutable vec3d vertex0 = inTriangle[0];
    immutable vec3d vertex1 = inTriangle[1];
    immutable vec3d vertex2 = inTriangle[2];

    vec3d edge1, edge2, h, s, q;

    double a, f, u, v;

    edge1 = vertex1 - vertex0;

    edge2 = vertex2 - vertex0;

    h = cross(rayVector, edge2);
    a = dot(edge1, h);

    if (a > -EPSILON && a < EPSILON)
        return false; // This ray is parallel to this triangle.

    f = 1.0 / a;

    s = rayOrigin - vertex0;

    u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    q = cross(s, edge1);

    v = f * dot(rayVector, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    immutable float t = f * dot(edge2, q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}


/************************************************************
 * Compute the associate Legendre polynomial
 * Params:
 *      l = degree
 *      m = order
 *      x = value 
 * Returns:
 *      double value of the Legendre Polynomial at x
 *
 * History: Adapted from Fortran 77 
 * (http://numerical.recipes/forum/showthread.php?t=200)
 *************************************************************/
pure nothrow double plegendre(in int l, in int m, in double x)
{
    int i, ll;
    double fact, oldfact, pll, pmm, pmmp1, omx2;

    pmm = 1.0;
    if (m > 0)
    {
        omx2 = (1.0 - x) * (1.0 + x);
        fact = 1.0;
        for (i = 1; i <= m; i++)
        {
            pmm *= omx2 * fact / (fact + 1.0);
            fact += 2.0;
        }
    }
    pmm = sqrt((2 * m + 1) * pmm / (4.0 * PI));
    if (m & 1)
        pmm = -pmm;
    if (l == m)
        return pmm;
    else
    {
        pmmp1 = x * sqrt(2.0 * m + 3.0) * pmm;
        if (l == (m + 1))
            return pmmp1;
        else
        {
            oldfact = sqrt(2.0 * m + 3.0);
            for (ll = m + 2; ll <= l; ll++)
            {
                fact = sqrt((4.0 * ll * ll - 1.0) / (ll * ll - m * m));
                pll = (x * pmmp1 - pmm / oldfact) * fact;
                oldfact = fact;
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}


/************************************************************
 * Compute the real Spherical Harmonics
 * Params:
 *      l = degree
 *      m = order
 *      theta = azimuthal angle
 *      phi = polar angle
 * Returns:
 *      double value of the Legendre Polynomial at x
 *************************************************************/
pure nothrow double ylm(in int l, in int m, in double theta, in double phi)
{
    const double m_sqrt2 = 1.41421356237309504880;

    double result;

    if (m < 0)
        result = m_sqrt2 * sin(-m * phi) * plegendre(l, -m, cos(theta));
    else if (m == 0)
        result = plegendre(l, m, cos(theta));
    else
        result = m_sqrt2 * cos(m * phi) * plegendre(l, m, cos(theta));

    return result;
}


/************************************************************
 * Compute the barycentric coordinates
 * Params:
 *      inTriangle = A, B & C vertices
 *      P = point for the coordinates
 * Returns:
 *      double[2] value for u & v
 * https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates 
 *************************************************************/
pure nothrow double[2] barycentric_coordinates(in vec3d[3] inTriangle, in vec3d P)
{
    immutable vec3d A = inTriangle[0];
    immutable vec3d B = inTriangle[1];
    immutable vec3d C = inTriangle[2];


    immutable double ABC = cross(B-A, C-A).length;

    immutable double CAP = cross(C-A, C-P).length;
    immutable double ABP = cross(A-B, A-P).length;

    immutable double u = CAP/ABC;
    immutable double v = ABP/ABC;

    return [u, v];
}
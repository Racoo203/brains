// Copyright (c) 2010 David W. Shattuck, PhD
// $Id: dfsexample.cpp 155 2010-08-24 22:15:30Z shattuck $
/**
 \brief  Example for use of the BrainSuite Surface Format file reader classes
 \author David Shattuck 
 \date   $Date: 2010-08-24 15:15:30 -0700 (Tue, 24 Aug 2010) $

 This example program takes one or two arguments. The first specifies an input filename.
 If only one argument is specified, then the program will read the input file and print
 out the first five elements of each vector present. If two arguments are provided, the
 program will read the input file and write it out to the location specified by the
 output file.
 */

#include "dfsurface.h"

//output function for UV coordinate class
inline std::ostream &operator<<(std::ostream &os, SILT::PointUV &t)
{
	return os<<t.u<<'\t'<<t.v;
}

//output function for Triangle  class
inline std::ostream &operator<<(std::ostream &os, SILT::Triangle &t)
{
	return os<<t.a<<'\t'<<t.b<<'\t'<<t.c;
}

//output function for 3D point class
inline std::ostream &operator<<(std::ostream &os, SILT::Point3D &p)
{
	return os<<p.x<<'\t'<<p.y<<'\t'<<p.z;
}

int main(int argc, char *argv[])
{
	SILT::DFSurface surface;

	// En este caso no se especifica el nombre del archivo de salida, 
    // por lo que se imprime en pantalla las indicaciones de uso.
	if (argc<2)
	{
		std::cerr<<"usage: dfssurface <dfs filename> { <dfs output file> }\n"<<std::endl;
		return 1;
	}
	std::string ifname(argv[1]);

	// Se lee el archivo de entrada. Si no se puede leer, 
    // se imprime un mensaje de error.
	if (surface.readDFS(ifname)==false)
	{
		std::cerr<<"error reading "<<ifname<<std::endl;
		return 1;
	}
	
	// Al corroborar que exista un archivo de entrada, se ejecuta el codigo de
    // la funcion main, que corresponde a lo que indica el paper.
	if (argc<3)
	{
		std::cout<<"read "<<ifname<<std::endl;
		std::cout<<(int)surface.triangles.size()<<" triangles."<<std::endl;
		std::cout<<(int)surface.vertices.size()<<" vertices."<<std::endl;
		std::cout<<(int)surface.vertexNormals.size()<<" vertex normals."<<std::endl;
		std::cout<<(int)surface.vertexColors.size()<<" vertex colors."<<std::endl;
		std::cout<<(int)surface.vertexUV.size()<<" vertex UVs."<<std::endl;
		std::cout<<(int)surface.vertexLabels.size()<<" vertex labels."<<std::endl;
		std::cout<<(int)surface.vertexAttributes.size()<<" vertex attributes."<<std::endl;
	
		const size_t limit = 5;
		const size_t nv = (surface.vertices.size()>limit) ? limit : surface.vertices.size();
		const size_t nt = (surface.triangles.size()>limit) ? limit : surface.triangles.size();
		std::cout<<"triangles\n";
		for (size_t i=0;i<nt;i++) std::cout<<surface.triangles[i]<<'\n';
		std::cout<<"vertices\n";
		for (size_t i=0;i<nt;i++) std::cout<<surface.vertices[i]<<'\n';
		if (surface.vertexColors.size()==surface.vertices.size())
		{
			std::cout<<"vertex colors\n";
			for (size_t i=0;i<nv;i++) std::cout<<surface.vertexColors[i]<<'\n';
		}
		if (surface.vertexNormals.size()==surface.vertices.size())
		{
			std::cout<<"vertex normals\n";
			for (size_t i=0;i<nv;i++) std::cout<<surface.vertexNormals[i]<<'\n';
		}
		if (surface.vertexUV.size()==surface.vertices.size())
		{
			std::cout<<"vertex UV\n";
			for (size_t i=0;i<nv;i++) std::cout<<surface.vertexUV[i]<<'\n';
		}
		if (surface.vertexLabels.size()==surface.vertices.size())
		{
			std::cout<<"vertex labels\n";
			for (size_t i=0;i<nv;i++) std::cout<<surface.vertexLabels[i]<<'\n';
		}
		if (surface.vertexAttributes.size()==surface.vertices.size())
		{
			std::cout<<"vertex attributes\n";
			for (size_t i=0;i<nv;i++) std::cout<<surface.vertexAttributes[i]<<'\n';
		}
	}

	// Se escribe el archivo de salida. Si no se puede escribir,
    // se imprime un mensaje de error.
	else
	{
		if (surface.writeDFS(argv[2])==false)
		{
			std::cerr<<"error writing "<<argv[2]<<std::endl;
			return 1;
		}
	}
	return 0;
}
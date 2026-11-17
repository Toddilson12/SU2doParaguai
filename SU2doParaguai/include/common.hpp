#pragma once
#include <cmath>

typedef unsigned long ulong;

namespace FVUtil
{
	enum class COORDINATES : ulong
	{
		X = 0,
		Y,
		Z
	};
	inline double sqr(double x)
	{
		return x * x;
	}
		
	inline double sign(double x)
	{
		return (x >= 0.0) ? 1.0 : 0.0;
	}


	struct Vec4
	{
		double v[4];
		
		Vec4 operator+(const Vec4& other) const
		{
			Vec4 result;
			for (size_t i = 0; i < 4; i++)
				result.v[i] = this->v[i] + other.v[i];
			return result;
		}

	
		Vec4 operator-(const Vec4& other) const
		{
			Vec4 result;
			for (size_t i = 0; i < 4; i++)
				result.v[i] = this->v[i] - other.v[i];
			return result;
		}

		Vec4 operator*(double scalar) const
		{
			Vec4 result;
			for (size_t i = 0; i < 4; i++)
				result.v[i] = this->v[i] * scalar;
			return result;
		}

		friend Vec4 operator*(double scalar, const Vec4& other)
		{
			return other * scalar;
		}

	};
	
	
	struct Vec2
	{
		double x, y;
		Vec2(double x_, double y_)
			: x(x_), y(y_)
		{
		}
		Vec2()
			: x(0.0), y(0.0)
		{
		}

		Vec2 operator-(const Vec2& other) const
		{
			return Vec2({ x - other.x, y - other.y });
		}
		
		Vec2 operator+(const Vec2& other) const
		{
			return Vec2({ x + other.x, y + other.y });
		}

		Vec2 operator*(double scalar) const
		{
			return Vec2({ x * scalar, y * scalar });
		}	

		
		static double GetDist(const Vec2& a, const Vec2& b)
		{
			return std::sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
		}

		static double Dot(const Vec2& a, const Vec2& b)
		{
			return a.x * b.x + a.y * a.y;
		}
		
		static double Mag(const Vec2& a)
		{
			return std::sqrt(sqr(a.x) + sqr(a.y));
		}

		static Vec2 Norm(const Vec2& a)
		{
			return a * (1.0 / Mag(a));
		}

		static double TriArea(const Vec2& a, const Vec2& b)
		{
			return 0.5 * std::abs(a.y * b.x - a.x * b.y);
		}
	};

	struct Block2DMatrix //contiguous matrix
	{
		double* block_ptr = nullptr;
		ulong rows = 0, cols = 0;
		
		Block2DMatrix operator=(double* other) //NOT MEMORY SAFE...
		{
			for (int n = 0; n < rows * cols; n++)
				block_ptr[n] = other[n];
			return *this;
		}
		
		Block2DMatrix operator=(const Block2DMatrix& other) //passes pointer of block_ptr basically
		{
			this->block_ptr = other.block_ptr;
			this->rows = other.rows;
			this->cols = other.cols;
			return *this;
		}

		double* operator()(ulong i)
		{
			return &block_ptr[i * cols];
		}

		double& operator()(ulong i, ulong j)
		{
			return block_ptr[i * cols + j];
		}

		static Block2DMatrix AllocateBlock2DMatrix(ulong rows, ulong cols)
		{
			Block2DMatrix mat;
			mat.rows = rows;
			mat.cols = cols;

			mat.block_ptr = new double[rows * cols];

			return mat;
		}

		void DeallocateMemory()
		{
			delete[] block_ptr;
			block_ptr = nullptr;
		}


	};
	
	

}
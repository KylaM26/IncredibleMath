#pragma once

#include <stdio.h>

namespace IncredibleMath {
	typedef struct Matrix2x2 {
		union {
			struct {
				float _11, _12,
					_21, _22;
			};

			float arr[4];
		};

		Matrix2x2 operator*(float scaler);

		inline float* operator[](const int index) {
			return &(arr[index * 2]);
		}

	} Matrix2x2, Mat2x2, M2x2;


	Matrix2x2 Matrix2x2::operator*(float scaler) {
		Matrix2x2 m;
		for (int i = 0; i < 2; i++)
			m.arr[i] = this->arr[i] * scaler;

		return m;
	}
}

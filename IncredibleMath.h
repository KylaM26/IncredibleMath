#ifndef INCREDIBLE_MATH
#define INCREDIBLE_MATH

#include <iostream>

namespace IncredibleMath {
#define PI 3.14159265
#define RAD_TO_DEG(f) (f * (180/PI)) 
#define DEG_TO_RAD(f) (f * (PI/180)) 

	/////////////////////////////////// VECTOR 2 //////////////////////////////
	typedef struct Vector2 {
		union {
			struct {
				float x, y;
			};

			float arr[2];
		};

		Vector2() {
			this->x = 0, this->y = 0;
		}

		Vector2(float x, float y) {
			this->x = x, this->y = y;
		}

		Vector2 operator+(const Vector2& v) {
			return Vector2(this->x + v.x, this->y + v.y);
		}

		Vector2 operator-(const Vector2& v) {
			return Vector2(this->x - v.x, this->y - v.y);
		}

		Vector2 operator*(const Vector2& v) {
			return Vector2(this->x * v.x, this->y * v.y);
		}

		Vector2 operator+(float f) {
			return Vector2(this->x + f, this->y + f);
		}

		Vector2 operator-(float f) {
			return Vector2(this->x - f, this->y - f);
		}

		Vector2 operator*(float f) {
			return Vector2(this->x * f, this->y * f);
		}

		float& operator[](const int index) {
			if (index < 2) return arr[index]; else printf("Index of: %d, is out of range!", index); return arr[0];
		}

		void operator=(Vector2 v) {
			this->x = v.x;
			this->y = v.y;
		}

		void operator+=(Vector2 v) {
			this->x += v.x, this->y += v.y;
		}

		void operator-=(Vector2 v) {
			this->x -= v.x, this->y -= v.y;
		}

		void operator*=(Vector2 v) {
			this->x *= v.x, this->y *= v.y;
		}

		void operator+=(float f) {
			this->x += f, this->y += f;
		}

		void operator-=(float f) {
			this->x -= f, this->y -= f;
		}

		void operator*=(float f) {
			this->x *= f, this->y *= f;
		}

		void operator++() {
			this->x++, this->y++;
		}

		void operator--() {
			this->x--; this->y--;
		}

		void Print() {
			printf("x: %f, y: %f\n", this->x, this->y);
		}
	} Vector2, Vec2, V2;

	//////////////////////////////// VECTOR 3 ///////////////////////////////
	typedef struct Vector3 {
		union {
			struct {
				float x, y, z;
			};

			struct {
				float r, g, b;
			};

			float arr[3];
		};

		Vector3() {
			this->x = 0, this->y = 0, this->z = 0;
		}

		Vector3(float x, float y, float z) {
			this->x = x, this->y = y, this->z = z;
		}

		Vector3 operator+(Vector3 v) {
			return Vector3(this->x + v.x, this->y + v.y, this->z + v.z);
		}

		Vector3 operator-(Vector3 v) {
			return Vector3(this->x - v.x, this->y - v.y, this->z - v.z);
		}

		Vector3 operator*(Vector3 v) {
			return Vector3(this->x * v.x, this->y * v.y, this->z * v.z);
		}

		Vector3 operator+(float f) {
			return Vector3(this->x + f, this->y + f, this->z + f);
		}

		Vector3 operator-(float f) {
			return Vector3(this->x - f, this->y - f, this->z - f);
		}

		Vector3 operator*(float f) {
			return Vector3(this->x * f, this->y * f, this->z * f);
		}

		float& operator[](const int index) {
			if (index < 3) return arr[index]; else printf("Index of: %d, is out of range!", index); return arr[0];
		}

		void operator=(Vector3 v) {
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
		}

		void operator+=(Vector3 v) {
			this->x += v.x, this->y += v.y, this->z += v.z;
		}

		void operator-=(Vector3 v) {
			this->x -= v.x, this->y -= v.y, this->z -= v.z;
		}

		void operator*=(Vector3 v) {
			this->x *= v.x, this->y *= v.y, this->z *= v.z;
		}

		void operator+=(float f) {
			this->x += f, this->y += f, this->z += f;
		}

		void operator-=(float f) {
			this->x -= f, this->y -= f, this->z -= f;
		}

		void operator*=(float f) {
			this->x *= f, this->y *= f, this->z *= f;
		}

		void operator++() {
			this->x++, this->y++, this->z++;
		}

		void operator--() {
			this->x--; this->y--, this->z--;
		}

		void Print() {
			printf("x: %f, y: %f, z: %f\n", this->x, this->y, this->z);
		}

	} Vector3, Vec3, V3;

	////////////////////////////// VECTOR 4 ///////////////////////////////
	typedef struct Vector4 {
		union {
			struct {
				float x, y, z, w;
			};

			struct {
				float r, g, b, a;
			};

			float arr[4];
		};

		Vector4() {
			this->x = 0, this->y = 0, this->z = 0;
		}
		Vector4(float x, float y, float z, float w) {
			this->x = x, this->y = y, this->z = z, this->w = w;
		}

		Vector4 operator+(const Vector4 v) {
			return Vector4(this->x + v.x, this->y + v.y, this->z + v.z, this->w + v.w);
		}

		Vector4 operator-(Vector4 v) {
			return Vector4(this->x - v.x, this->y - v.y, this->z - v.z, this->w - v.w);
		}

		Vector4 operator*(Vector4 v) {
			return Vector4(this->x * v.x, this->y * v.y, this->z * v.z, this->w * v.w);
		}

		Vector4 operator+(float f) {
			return Vector4(this->x + f, this->y + f, this->z + f, this->w + f);
		}

		Vector4 operator-(float f) {
			return Vector4(this->x - f, this->y - f, this->z - f, this->w - f);
		}

		Vector4 operator*(float f) {
			return Vector4(this->x * f, this->y * f, this->z * f, this->w * f);
		}

		float& operator[](const int index) {
			if (index < 4) return arr[index]; else printf("Index of: %d, is out of range!", index); return arr[0];
		}

		void operator=(Vector4 v) {
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
			this->w = v.w;
		}

		void operator+=(Vector4 v) {
			this->x += v.x, this->y += v.y, this->z += v.z, this->w += v.w;
		}

		void operator-=(Vector4 v) {
			this->x -= v.x, this->y -= v.y, this->z -= v.z, this->w -= v.w;
		}

		void operator*=(Vector4 v) {
			this->x *= v.x, this->y *= v.y, this->z *= v.z, this->w *= v.w;
		}

		void operator+=(float f) {
			this->x += f, this->y += f, this->z += f, this->w += f;
		}

		void operator-=(float f) {
			this->x -= f, this->y -= f, this->z -= f, this->w -= f;
		}

		void operator*=(float f) {
			this->x *= f, this->y *= f, this->z *= f, this->w *= f;
		}

		void operator++() {
			this->x++, this->y++, this->z++, this->w++;
		}

		void operator--() {
			this->x--; this->y--, this->z--, this->w--;
		}

		void Print() {
			printf("x: %f, y: %f, z: %f, w: %f\n", this->x, this->y, this->z, this->w);
		}

	} Vector4, Vec4, V4;

	///////////////////////// MATRIX2x2 //////////////////////////////
	typedef struct Matrix2x2 {
		union {
			struct {
				float _11, _12,
					_21, _22;
			};

			float arr[4];
		};

		Matrix2x2() {
			this->_11 = 0.f, this->_12 = 0.f, this->_21 = 0.f, this->_22 = 0.f;
		}


		Matrix2x2(float _11, float _12, float _21, float _22) {
			this->_11 = _11, this->_12 = _12, this->_21 = _21, this->_22 = _22;
		}

		Matrix2x2 operator*(float scaler) {
			Matrix2x2 m;
			for (int i = 0; i < 2; i++)
				m.arr[i] = this->arr[i] * scaler;

			return m;
		}

		inline float* operator[](const int index) {
			return &(arr[index * 2]);
		}

	} Matrix2x2, Mat2x2, M2x2;

	///////////////////////// MATRIX3x3 //////////////////////////////
	typedef struct Matrix3x3 {
		union {
			struct {
				float _11, _12, _13,
					_21, _22, _23,
					_31, _32, _33;
			};

			float arr[9];
		};

		Matrix3x3 operator*(float scaler) {
			Matrix3x3 m;
			for (int i = 0; i < 3; i++)
				m.arr[i] = this->arr[i] * scaler;

			return m;
		}

		inline float* operator[](const int index) {
			return &(arr[index * 3]);
		}

	} Matrix3x3, Mat3x3, M3x3;

	///////////////////////// MATRIX34x4 //////////////////////////////
	typedef struct Matrix4x4 {
		union {
			struct {
				float _11, _12, _13, _14,
					_21, _22, _23, _24,
					_31, _32, _33, _34,
					_41, _42, _43, _44;
			};

			float arr[16];
		};

		Matrix4x4() {
			this->_11 = 1.0f; this->_12 = 0.0f; this->_13 = 0.0f; this->_14 = 0.0f;
			this->_21 = 0.0f; this->_22 = 1.0f; this->_23 = 0.0f; this->_24 = 0.0f;
			this->_31 = 0.0f; this->_32 = 0.0f; this->_33 = 1.0f; this->_34 = 0.0f;
			this->_11 = 0.0f; this->_42 = 0.0f; this->_43 = 0.0f; this->_44 = 1.0f;
		}

		Matrix4x4(float xPos, float yPos, float zPos) {
			this->_11 = 1.0f; this->_12 = 0.0f; this->_13 = 0.0f; this->_14 = xPos;
			this->_21 = 0.0f; this->_22 = 1.0f; this->_23 = 0.0f; this->_24 = yPos;
			this->_31 = 0.0f; this->_32 = 0.0f; this->_33 = 1.0f; this->_34 = zPos;
			this->_11 = 0.0f; this->_42 = 0.0f; this->_43 = 0.0f; this->_44 = 1.0f;
		}

		Matrix4x4(float _11, float _12, float _13, float _14,
			float _21, float _22, float _23, float _24,
			float _31, float _32, float _33, float _34,
			float _41, float _42, float _43, float _44) {

			this->_11 = _11; this->_12 = _12; this->_13 = _13; this->_14 = _14;
			this->_21 = _21; this->_22 = _22; this->_23 = _23; this->_24 = _24;
			this->_31 = _34; this->_32 = _32; this->_33 = _33; this->_34 = _34;
			this->_11 = _41; this->_42 = _42; this->_43 = _43; this->_44 = _44;
		}


		Matrix4x4 operator*(Matrix4x4& matrix) {
			Matrix4x4 m;

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					m[i][j] =
						*this[i][0] * matrix[0][j] +
						*this[i][1] * matrix[1][j] +
						*this[i][2] * matrix[2][j] +
						*this[i][3] * matrix[3][j];
				}
			}

			return m;
		}

		Matrix4x4 operator*(float scaler) {
			Matrix4x4 m;
			for (int i = 0; i < 4; i++)
				m.arr[i] = this->arr[i] * scaler;

			return m;
		}

		inline float* operator[](const int index) {
			return &(arr[index * 4]);
		}

	} Matrix4x4, Mat4x4, M4x4;

	/////////////////////// MATH OPERATIONS ///////////////////////////
		// DOT
	inline float Dot(Vector2 a, Vector2 b) { // This function does not normalize the vector.
		return ((a.x * b.x) + (a.y * b.y));
	}

	inline float Dot(Vector3 a, Vector3 b) { // This function does not normalize the vector.
		return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
	}

	inline float Dot(Vector4 a, Vector4 b) {
		return ((a.x * b.x) + (a.y * b.y));
	}

	// MAGNITUTDE
	inline float Magnitude(Vector2 v) {
		return sqrtf(Dot(v, v));
	}

	inline float Magnitude(Vector3 v) {
		return sqrtf(Dot(v, v));
	}

	inline float Magnitude(Vector4 v) {
		return sqrtf(Dot(v, v));
	}

	inline float MagnitudeSq(Vector2 v) { // Instead of using the normal Magnitutde function just used this and multiply the value it is being checked against twice. To avoid using the sqrt() function
		return Dot(v, v);
	}

	inline float MagnitudeSq(Vector3 v) { // Instead of using the normal Magnitutde function just used this and multiply the value it is being checked against twice. To avoid using the sqrt() function
		return Dot(v, v);
	}

	inline float MagnitudeSq(Vector4 v) {
		return Dot(v, v);
	}

	inline Vector2 Normalize(Vector2 v) {
		return (v * (1.f / Magnitude(v)));
	}

	inline Vector3 Normalize(Vector3 v) {
		return (v * (1.f / Magnitude(v)));
	}

	inline Vector4 Normalize(Vector4 v) {
		return (v * (1.f / Magnitude(v)));
	}

	// Angle
	inline float GetAngleBetweenVectors(Vector2 a, Vector2 b) {
		return acosf(Dot(a, b));
	}

	inline float GetAngleBetweenVectors(Vector3 a, Vector3 b) {
		return acosf(Dot(a, b));
	}

	inline float GetAngleBetweenVectors(Vector4 a, Vector4 b) {
		return (acosf(Dot(a, b)));
	}

	// CROSS
	inline Vector3 Cross(const Vector3 a, const Vector3 b) {
		Vector3 result = Vector3();
		result.x = a.y * b.z - a.z * b.y;
		result.y = a.z * b.x - a.x * b.z;
		result.z = a.x * b.y - a.y * b.x;
		return result;
	}

	inline Vector2 Project(Vector2 length, Vector2 direction) { // B will be projected onto A
		float dP = Dot(length, direction);
		float magSq = MagnitudeSq(direction);
		return direction * (dP / magSq);
	}

	inline Vector3 Project(Vector3 length, Vector3 direction) { // B will be projected onto A
		float dP = Dot(length, direction);
		float magSq = MagnitudeSq(direction);
		return direction * (dP / magSq);
	}

	inline Vector4 Project(Vector4 length, Vector4 direction) {
		float dP = Dot(length, direction);
		float magSq = MagnitudeSq(direction);
		return direction * (dP / magSq);
	}

	inline Vector2 Perpendicular(Vector2 length, Vector2 direction) {
		return length - Project(length, direction);
	}

	inline Vector3 Perpendicular(Vector3 length, Vector3 direction) {
		return length - Project(length, direction);
	}

	inline Vector4 Perpendicular(Vector4& length, Vector4& direction) {
		return length - Project(length, direction);
	}

	// REFLECTION
	inline Vector2 Reflect(Vector2 v, Vector2 normal) {
		float dot = Dot(v, normal);
		return (v - normal * (dot * 2.f));
	}

	inline Vector3 Reflect(Vector3 v, Vector3 normal) {
		float dot = Dot(v, normal);
		return (v - normal * (dot * 2.f));
	}

	inline Vector4 Reflect(Vector4& v, Vector4& normal) {
		float dot = Dot(v, normal);
		return (v - normal * (dot * 2.f));
	}

	// MATRIX OPERATRIONS
	inline Mat4x4 LookAt(Vector3 position, Vector3 target, Vector3 up) {
		Vector3 forward = Normalize(target - position);
		Vector3 right = Normalize(Cross(up, forward));
		Vector3 newUp = Cross(forward, right);

		return Mat4x4(
			right.x, newUp.x, forward.x, 0.0f,
			right.y, newUp.y, forward.y, 0.0f,
			right.z, newUp.z, forward.z, 0.0f,
			-Dot(right, position),
			-Dot(newUp, position),
			-Dot(forward, position), 1.0f
		);
	}

	inline Mat4x4 Projection(float fov, float aspect, float zNear, float zFar) {
		float tanHalfFov = tanf(DEG_TO_RAD(fov * 0.5f));
		float fovY = 1.f / tanHalfFov;
		float fovX = fovY / aspect;

		Mat4x4 result;
		result._11 = fovX;
		result._22 = fovY;
		result._33 = zFar / (zFar - zNear);
		result._34 = 1.0f;
		result._43 = -zNear * result._33;
		result._44 = 0.0f;

		return result;
	}

	inline Mat4x4 Ortho(float left, float right, float bottom, float top, float zNear, float zFar) {
		float _11 = 2.0f / (right - left);
		float _22 = 2.0f / (top - bottom);
		float _33 = 1.0f / (zFar - zNear);
		float _41 = (left + right) / (left - right);
		float _42 = (top + bottom) / (bottom - top);
		float _43 = zFar / (zNear - zFar);

		return Mat4x4(
			_11, 0.0f, 0.0f, 0.0f,
			0.0f, _22, 0.0f, 0.0f,
			0.0f, 0.0f, _33, 0.0f,
			_41, _42, _43, 1.0f
		);
	}

	//inline Mat4x4 Transform(Vector3 scale, Vector3 eulerRot, Vector3 translate) {
	//	return Scale()
	//}
}

#endif // !
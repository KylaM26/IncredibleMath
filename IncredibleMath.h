#pragma once
#include <iostream>

#define PI 3.14159265
#define RAD_TO_DEG(f) (f * (180/PI)) 
#define DEG_TO_RAD(f) (f * (PI/180)) 

#define CMP(x, y) (fabsf((x) – (y)) <= FLT_EPSILON * fmaxf(1.0f, fmaxf(fabsf(x), fabsf(y))))

namespace IncredibleMath {
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


	typedef struct Matrix4x4 {
		union {
			float m[4][4];
			float a[16];
		};

		Matrix4x4() { // IDENTITY
			m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = 0.0f;
			m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = 0.0f;
			m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = 0.0f;
			m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
		}

		Matrix4x4(float a, float b, float c) { // POSITION
			m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = a;
			m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = b;
			m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = c;
			m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
		}

		Matrix4x4(
			float _11, float _12, float _13, float _14,
			float _21, float _22, float _23, float _24,
			float _31, float _32, float _33, float _34,
			float _41, float _42, float _43, float _44
		) {
			m[0][0] = _11; m[0][1] = _12; m[0][2] = _13; m[0][3] = _14;
			m[1][0] = _21; m[1][1] = _22; m[1][2] = _23; m[1][3] = _24;
			m[2][0] = _31; m[2][1] = _32; m[2][2] = _33; m[2][3] = _34;
			m[3][0] = _41; m[3][1] = _42; m[3][2] = _43; m[3][3] = _44;

		}

		inline Matrix4x4 operator=(const Matrix4x4& e) {
			this->m[0][0] = e.m[0][0]; this->m[0][1] = e.m[0][1]; this->m[0][2] = e.m[0][2]; this->m[0][3] = e.m[0][3];
			this->m[1][0] = e.m[1][0]; this->m[1][1] = e.m[1][1]; this->m[1][2] = e.m[1][2]; this->m[1][3] = e.m[1][3];
			this->m[2][0] = e.m[2][0]; this->m[2][1] = e.m[2][1]; this->m[2][2] = e.m[2][2]; this->m[2][3] = e.m[2][3];
			this->m[3][0] = e.m[3][0]; this->m[3][1] = e.m[3][1]; this->m[3][2] = e.m[3][2]; this->m[3][3] = e.m[3][3];
			return *this;
		}

		inline Matrix4x4 operator*(const Matrix4x4& matrix) {
			Mat4x4 f;
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					f.m[i][j] =
						this->m[i][0] * matrix.m[0][j] +
						this->m[i][1] * matrix.m[1][j] +
						this->m[i][2] * matrix.m[2][j] +
						this->m[i][3] * matrix.m[3][j];
				}
			} return f;
		}

		inline Matrix4x4 operator*(float scaler) {
			for (int i = 0; i < 4; i++)
				this->a[i] = this->a[i] * scaler;
			return *this;
		}

		inline float* operator[](const int i) {
			return &(a[i * 4]);
		}

		inline void Translate(float x, float y, float z) {
			m[0][0] = 0.0f;	m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = x;
			m[1][0] = 0.0f; m[1][1] = 0.0f;	m[1][2] = 0.0f; m[1][3] = y;
			m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 0.0f;	m[2][3] = z;
			m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
		}

		inline void Translate(const Vector3& v) {
			m[0][0] = 0.0f;	m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = v.x;
			m[1][0] = 0.0f; m[1][1] = 0.0f;	m[1][2] = 0.0f; m[1][3] = v.y;
			m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 0.0f;	m[2][3] = v.z;
			m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
		}

	} M4x4, Mat4x4;

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
		result.m[0][0] = fovX;
		result.m[1][1] = fovY;
		result.m[2][2] = zFar / (zFar - zNear);
		result.m[2][3] = 1.0f;
		result.m[3][2] = -zNear * result.m[2][2];
		result.m[3][3] = 0.0f;

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

	inline void Transpose(const float* srcMat, float* dstMat, int srcRows, int srcCols) {
		for (int i = 0; i < srcRows * srcCols; i++) {
			int row = i / srcRows;
			int col = i % srcRows;
			dstMat[i] = srcMat[srcCols * col + row];
		}
	}

	inline Mat4x4 Transpose(const Mat4x4& matrix) {
		Mat4x4 result;
		Transpose(matrix.a, result.a, 4, 4);
		return result;
	}

	inline Mat4x4 Scale(float x, float y, float z) {
		return Mat4x4(
			x, 0.0f, 0.0f, 0.0f,
			0.0f, y, 0.0f, 0.0f,
			0.0f, 0.0f, z, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f);
	}


	inline Mat4x4 XRotation(float angle) {
		angle = DEG_TO_RAD(angle);
		return Mat4x4(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, cosf(angle), -sinf(angle), 0.0f,
			0.0f, sinf(angle), cosf(angle), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f);
	}

	inline Mat4x4 YRotation(float angle) {
		angle = DEG_TO_RAD(angle);
		return Mat4x4(
			cosf(angle), 0.0f, sinf(angle), 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			-sinf(angle), 0.0f, cosf(angle), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f);
	}

	inline Mat4x4 ZRotation(float angle) {
		angle = DEG_TO_RAD(angle);
		return Mat4x4(
			cosf(angle), -sinf(angle), 0.0f, 0.0f,
			sinf(angle), cosf(angle), 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f);
	}

	inline Mat4x4 Rotate(float yaw, float pitch, float roll) {
		return ZRotation(roll) * XRotation(pitch) * YRotation(yaw);
	}


}
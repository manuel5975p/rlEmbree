#ifndef RAYMATH_OPS
#define RAYMATH_OPS
#include "raylib.h"

//------------------Vector2-----------------//
inline Vector2 operator+ (Vector2 lhs, const Vector2& rhs)
{
	return Vector2{ lhs.x + rhs.x, lhs.y + rhs.y };
}

inline Vector2 operator+ (Vector2 lhs, const float& rhs)
{
	return Vector2{ lhs.x + rhs, lhs.y + rhs };
}

inline Vector2 operator- (Vector2 lhs, const Vector2& rhs)
{
	return Vector2{ lhs.x - rhs.x, lhs.y - rhs.y };
}

inline Vector2 operator- (Vector2 lhs, const float& rhs)
{
	return Vector2{ lhs.x - rhs, lhs.y - rhs };
}

inline Vector2 operator* (Vector2 lhs, const float& rhs)
{
	return Vector2{ lhs.x  * rhs, lhs.y * rhs };
}

inline Vector2 operator* (Vector2 lhs, const Vector2& rhs)
{
	return Vector2{ lhs.x * rhs.x, lhs.y * rhs.y };
}

inline Vector2 operator/ (Vector2 lhs, const float& rhs)
{
	return Vector2{ lhs.x / rhs, lhs.y / rhs };
}

inline Vector2 operator/ (Vector2 lhs, const Vector2& rhs)
{
	return Vector2{ lhs.x / rhs.y, lhs.y / rhs.y };
}

//------------------Vector3-----------------//
inline Vector3 operator+ (Vector3 lhs, const Vector3& rhs)
{
	return Vector3{ lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
}

inline Vector3 operator+ (Vector3 lhs, const float& rhs)
{
	return Vector3{ lhs.x + rhs, lhs.y + rhs, lhs.z + rhs };
}

inline Vector3 operator- (Vector3 lhs, const Vector3& rhs)
{
	return Vector3{ lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
}

inline Vector3 operator- (Vector3 lhs, const float& rhs)
{
	return Vector3{ lhs.x - rhs, lhs.y - rhs, lhs.z - rhs };
}

inline Vector3 operator* (Vector3 lhs, const float& rhs)
{
	return Vector3{ lhs.x * rhs, lhs.y * rhs, lhs.z * rhs };
}

inline Vector3 operator* (Vector3 lhs, const Vector3& rhs)
{
	return Vector3{ lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z };
}

inline Vector3 operator/ (Vector3 lhs, const float& rhs)
{
	return Vector3{ lhs.x / rhs, lhs.y / rhs, lhs.z / rhs };
}

inline Vector3 operator/ (Vector3 lhs, const Vector3& rhs)
{
	return Vector3{ lhs.x / rhs.y, lhs.y / rhs.y, lhs.z / rhs.z };
}

static constexpr Vector3 Vector3X = { 1, 0, 0 };
static constexpr Vector3 Vector3Y = { 0, 1, 0 };
static constexpr Vector3 Vector3Z = { 0, 0, 1 };
template<typename str>
inline str& operator<<(str& s, Vector2 vec){
	return s << "{" << vec.x << ", " << vec.y << "}";
}
template<typename str>
inline str& operator<<(str& s, Vector3 vec){
	return s << "{" << vec.x << ", " << vec.y << ", " << vec.z << "}";
}

#endif
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[derive(Clone, Copy, Default, PartialEq)]
pub struct Vec2 {
    pub x: f32,
    pub y: f32,
}

impl Vec2 {
    pub const ONE: Vec2 = Vec2::new(1.0, 1.0);
    pub const ZERO: Vec2 = Vec2::new(0.0, 0.0);

    pub const fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn dot(&self, rhs: Vec2) -> f32 {
        self.x * rhs.x + self.y * rhs.y
    }

    pub fn cross(&self, rhs: Vec2) -> f32 {
        Mat2::new(*self, rhs).det()
    }

    pub fn cross_scalar(lhs: f32, rhs: Vec2) -> Vec2 {
        Vec2::new(-lhs * rhs.y, lhs * rhs.x)
    }

    pub fn cross_scalar2(lhs: Vec2, rhs: f32) -> Vec2 {
        Vec2::new(lhs.x * rhs, lhs.y * -rhs)
    }

    pub(crate) fn length_squared(&self) -> f32 {
        self.dot(*self)
    }

    pub(crate) fn length(&self) -> f32 {
        self.length_squared().sqrt()
    }

    pub fn normalized(&self) -> Self {
        *self / self.length()
    }
}

impl Debug for Vec2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("Vec2").field(&self.x).field(&self.y).finish()
    }
}

impl Neg for Vec2 {
    type Output = Vec2;
    fn neg(self) -> Self::Output {
        Vec2::new(-self.x, -self.y)
    }
}

impl Add<Vec2> for Vec2 {
    type Output = Vec2;
    fn add(self, rhs: Vec2) -> Self::Output {
        Vec2::new(self.x + rhs.x, self.y + rhs.y)
    }
}
impl AddAssign<Vec2> for Vec2 {
    fn add_assign(&mut self, rhs: Vec2) {
        *self = *self + rhs;
    }
}

impl Sub<Vec2> for Vec2 {
    type Output = Vec2;
    fn sub(self, rhs: Vec2) -> Self::Output {
        Vec2::new(self.x - rhs.x, self.y - rhs.y)
    }
}
impl SubAssign<Vec2> for Vec2 {
    fn sub_assign(&mut self, rhs: Vec2) {
        *self = *self - rhs;
    }
}

impl Mul<f32> for Vec2 {
    type Output = Vec2;
    fn mul(self, rhs: f32) -> Self::Output {
        Vec2::new(self.x * rhs, self.y * rhs)
    }
}
impl MulAssign<f32> for Vec2 {
    fn mul_assign(&mut self, rhs: f32) {
        *self = *self * rhs;
    }
}

impl Mul<Vec2> for f32 {
    type Output = Vec2;
    fn mul(self, rhs: Vec2) -> Self::Output {
        Vec2::new(self * rhs.x, self * rhs.y)
    }
}

impl Div<f32> for Vec2 {
    type Output = Vec2;
    fn div(self, rhs: f32) -> Self::Output {
        Vec2::new(self.x / rhs, self.y / rhs)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Mat2 {
    pub m0_col: Vec2,
    pub m1_col: Vec2,
}

impl Mat2 {
    pub const IDENTITY: Mat2 = Mat2 {
        m0_col: Vec2::new(1., 0.),
        m1_col: Vec2::new(0., 1.),
    };
    pub fn new(m0: Vec2, m1: Vec2) -> Self {
        Self {
            m0_col: m0,
            m1_col: m1,
        }
    }

    pub fn rotation(angle: f32) -> Self {
        Self::from([angle.cos(), -angle.sin(), angle.sin(), angle.cos()])
    }

    pub fn to_array(&self) -> [f32; 4] {
        [self.m0_col.x, self.m0_col.y, self.m1_col.x, self.m1_col.y]
    }

    pub fn det(&self) -> f32 {
        self.m0_col.x * self.m1_col.y - self.m0_col.y * self.m1_col.x
    }
}

impl From<[f32; 4]> for Mat2 {
    fn from(value: [f32; 4]) -> Self {
        Self::new(Vec2::new(value[0], value[1]), Vec2::new(value[2], value[3]))
    }
}

impl Mul<Vec2> for Mat2 {
    type Output = Vec2;

    fn mul(self, rhs: Vec2) -> Self::Output {
        Vec2::new(
            self.m0_col.x * rhs.x + self.m0_col.y * rhs.y,
            self.m1_col.x * rhs.x + self.m1_col.y * rhs.y,
        )
    }
}

// #[derive(Clone, Copy, Debug)]
// pub struct Mat3 {
//     m00: f32,
//     m01: f32,
//     m02: f32,
//     m10: f32,
//     m11: f32,
//     m12: f32,
//     m20: f32,
//     m21: f32,
//     m22: f32,
// }

// impl Mat3 {
//     pub const IDENTITY: Mat3 = Mat3 {
//         m00: 1.,
//         m01: 0.,
//         m02: 0.,
//         m10: 0.,
//         m11: 1.,
//         m12: 0.,
//         m20: 0.,
//         m21: 0.,
//         m22: 1.,
//     };
//     pub fn new(
//         m00: f32,
//         m01: f32,
//         m02: f32,
//         m10: f32,
//         m11: f32,
//         m12: f32,
//         m20: f32,
//         m21: f32,
//         m22: f32,
//     ) -> Self {
//         Self {
//             m00,
//             m01,
//             m02,
//             m10,
//             m11,
//             m12,
//             m20,
//             m21,
//             m22,
//         }
//     }

//     pub fn translation_rotation(translate: Vec2, rotate: f32) -> Self {
//         Mat3::new(rotate.cos(), -rotate.sin(), 0., rotate.sin(), rotate.cos(), 0., 0., 0., 1.) * Mat3::new(1., 0., translate.x,0.,1.,translate.y,)

//         Self {
//             m00: rotate.cos(),
//             m01: -rotate.sin(),
//             m02: translate.x * rotate.cos() - translate.y * rotate.sin(),
//             m10: rotate.sin(),
//             m11: rotate.cos(),
//             m12: translate.x * rotate.sin() + translate.y * rotate.cos(),
//             m20: 0.,
//             m21: 0.,
//             m22: 1.,
//         }
//     }

//     pub fn determinant(&self) -> f32 {
//         self.m00 * (self.m11 * self.m22 - self.m12 * self.m21)
//             - self.m01 * (self.m10 * self.m22 - self.m12 * self.m20)
//             + self.m02 * (self.m10 * self.m21 - self.m11 * self.m20)
//     }

//     pub fn inverse(&self) -> Option<Mat3> {
//         let det = self.m00 * (self.m11 * self.m22 - self.m12 * self.m21)
//             - self.m01 * (self.m10 * self.m22 - self.m12 * self.m20)
//             + self.m02 * (self.m10 * self.m21 - self.m11 * self.m20);

//         if det.abs() < 1e-6 {
//             return None;
//         }

//         let inv_det = 1.0 / det;

//         Some(Mat3 {
//             m00: (self.m11 * self.m22 - self.m12 * self.m21) * inv_det,
//             m01: (self.m02 * self.m21 - self.m01 * self.m22) * inv_det,
//             m02: (self.m01 * self.m12 - self.m02 * self.m11) * inv_det,
//             m10: (self.m12 * self.m20 - self.m10 * self.m22) * inv_det,
//             m11: (self.m00 * self.m22 - self.m02 * self.m20) * inv_det,
//             m12: (self.m02 * self.m10 - self.m00 * self.m12) * inv_det,
//             m20: (self.m10 * self.m21 - self.m11 * self.m20) * inv_det,
//             m21: (self.m01 * self.m20 - self.m00 * self.m21) * inv_det,
//             m22: (self.m00 * self.m11 - self.m01 * self.m10) * inv_det,
//         })
//     }

//     pub(crate) fn no_translation(&self) -> Self {
//         Self {
//             m02: 0.,
//             m12: 0.,
//             ..self.clone()
//         }
//     }
// }

// impl Mul<Mat3> for Mat3 {
//     type Output = Mat3;

//     fn mul(self, rhs: Mat3) -> Self::Output {
//         Mat3 {
//             m00: self.m00 * rhs.m00 + self.m01 * rhs.m10 + self.m02 * rhs.m20,
//             m01: self.m00 * rhs.m01 + self.m01 * rhs.m11 + self.m02 * rhs.m21,
//             m02: self.m00 * rhs.m02 + self.m01 * rhs.m12 + self.m02 * rhs.m22,
//             m10: self.m10 * rhs.m00 + self.m11 * rhs.m10 + self.m12 * rhs.m20,
//             m11: self.m10 * rhs.m01 + self.m11 * rhs.m11 + self.m12 * rhs.m21,
//             m12: self.m10 * rhs.m02 + self.m11 * rhs.m12 + self.m12 * rhs.m22,
//             m20: self.m20 * rhs.m00 + self.m21 * rhs.m10 + self.m22 * rhs.m20,
//             m21: self.m20 * rhs.m01 + self.m21 * rhs.m11 + self.m22 * rhs.m21,
//             m22: self.m20 * rhs.m02 + self.m21 * rhs.m12 + self.m22 * rhs.m22,
//         }
//     }
// }

// impl Mul<Vec2> for Mat3 {
//     type Output = Vec2;

//     fn mul(self, rhs: Vec2) -> Self::Output {
//         Vec2::new(
//             self.m00 * rhs.x + self.m01 * rhs.y + self.m02,
//             self.m10 * rhs.x + self.m11 * rhs.y + self.m12,
//         )
//     }
// }

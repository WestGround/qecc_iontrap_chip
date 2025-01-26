pub mod circuit;
pub mod schedule;
pub mod error_generator;

pub const SQ_TIME: usize = 1;
pub const TQ_TIME: usize = 5;
pub const SWAP_TIME: usize = 3*TQ_TIME;
pub const SHUTTLE_TIME: usize = 15;
pub const MEASURE_TIME: usize = 1;

pub const C7_QECTIME: usize = 5*TQ_TIME;
pub const C7_CODE: (usize, usize, usize) = (7, 1, 3);
pub const C17_QECTIME: usize = 10*TQ_TIME;
pub const C17_CODE: (usize, usize, usize) = (17, 1, 5);
pub const C31_QECTIME: usize = 10*TQ_TIME;
pub const C31_CODE: (usize, usize, usize) = (31, 1, 7);
pub const CODE_TIMES: [usize; 3] = [C7_QECTIME, C17_QECTIME, C31_QECTIME];
pub const CODES: [(usize, usize, usize); 3] = [C7_CODE, C17_CODE, C31_CODE];

pub const SQ_ERR: f64 = 0.1;
pub const TQ_ERR: f64 = 1.0;
pub const SWAP_ERR: f64 = 3.0;
pub const SHUTTLE_ERR: f64 = 0.1;
pub const MEASURE_ERR: f64 = 0.1;
pub const QEC_ERR: f64 = 5.0;
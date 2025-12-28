// src/rand.rs
//========================================================================================//
//                                                                                        //
// NAME  : rand.rs                                                                        //
// VER.  : 0.0.0.1                                                                        //
//                                                                                        //
// TYPE  : Rust Source File                                                               //
//                                                                                        //
// AUTHOR: Dylan Jurski                                                                   //
//                                                                                        //
// DESCRIPTION: Random number generation derived from pcg-random for C.                   //
// Original code is available at https://github.com/imneme/pcg-c-basic/.                  //
// For more infromation, visit https://pcg-random.org/.                                   //
//                                                                                        //
//========================================================================================//
//                                                                                        //
// Input parameters.                                                                      //
//                                                                                        //
//========================================================================================//
//                                                                                        //
// None                                                                                   //
//                                                                                        //
//========================================================================================//
//                                                                                        //
// Use Rust crates.                                                                       //
//                                                                                        //
//========================================================================================//
//                                                                                        //
// Define constants.                                                                      //
//                                                                                        //
//========================================================================================//
//                                                                                        //
// Define global variables.                                                               //
//                                                                                        //
//========================================================================================//
pub struct PCG32 {
    state: u64, // RNG state.  All values are possible.
    inc: u64,   // Controls which RNG sequence (stream) is
                // selected. Must *always* be odd.
}

pub type Prn32 = PCG32;
//========================================================================================//
//                                                                                        //
// Helper functions and impls.                                                            //
//                                                                                        //
//========================================================================================//
impl Prn32 {
    pub fn new(thread_id: u64) -> Self {
        let seed = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos() as u64;

        // Initialize in zero state with stream based on thread ID
        let mut rng = Self {
            state: 0,
            inc: (thread_id << 1) | 1,
        };

        // Move away from 0
        rng.next_u32();

        // Perform warmup to mix the seed into the state
        rng.state = seed.wrapping_add(rng.inc);

        // Advance the state once to avoid initial correlation and return u32
        rng.next_u32();

        rng
    }

    #[inline]
    fn next_u32(&mut self) -> u32 {
        let old_state = self.state;

        // Linear Congruential Generator step
        // We use the wrapping_mul to allow the 64-bit overflow naturally
        self.state = old_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(self.inc);

        let xorshifted = (((old_state >> 18) ^ old_state) >> 27) as u32;
        xorshifted.rotate_right((old_state >> 59) as u32)
    }
    /*
    #[inline]
    pub fn next_bitshift_f32(&mut self) -> f32 {
        // The bit-hack: range [0.0, 1.0)
        let u = self.next_u32();
        f32::from_bits((u >> 9) | 0x3f800000) - 1.0
    }
     */
    #[inline]
    pub fn next_f32(&mut self) -> f32 {
        // Multiply by 1.0 / 2^32 or 2.3283064E-10
        // This ensures the result is in range [0.0, 1.0)
        let u = self.next_u32();
        (u as f32) * (2.3283064E-10)
    }

    #[inline]
    pub fn next_f32_range(&mut self, min: f32, max: f32) -> f32 {
        let range = max - min;
        let rand_f32 = self.next_f32();
        rand_f32.mul_add(range, min)
    }
}

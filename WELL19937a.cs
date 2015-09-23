/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

using System;
using System.Runtime.InteropServices;

namespace RandomAlgorithms.WELL
{
	[Serializable]
	[ComVisible(true)]
	public class WELL19937a : Random
	{
		private static double FACT = 2.32830643653869628906e-10;
		private uint	m_i;
		private uint[]	m_state = new uint[624];
		private uint	m_z0;
		private uint	m_z1;
		private uint	m_z2;
		private uint	m_y;
		private int		m_case = 1;
		
		public WELL19937a() : this(DateTime.Now.Millisecond)
		{ }
		public WELL19937a(int seed)
		{
			var rnd = new Random(seed);
			for (int i = 0; i < 624; ++i)
				 m_state[i] = unchecked((uint)rnd.Next());
		}
		public WELL19937a(uint[] seed)
		{
			if (seed == null)
				throw new ArgumentNullException("The seed must not be null.");
			if (seed.Length != 624)
				throw new ArgumentException("Length of seed must be 624.");

			for (int i = 0; i < 624; ++i)
				 m_state[i] = seed[i];
		}

		private static uint MAT0POS(int t, uint v)
		{
			return v ^ (v >> t);
		}
		private static uint MAT0NEG(int t, uint v)
		{
			return v ^ (v << -t);
		}
		private static uint MAT3POS(int t, uint v)
		{
			return v >> t;
		}

		/// To obtain the WELL19937c
		public bool TEMPERING { get; set; }

		protected override double Sample()
		{
			switch (m_case)
			{
				case 1: return Case1();
				case 2: return Case2();
				case 3: return Case3();
				case 4: return Case4();
				case 5: return Case5();
				case 6: return Case6();
			}

			return 0;
		}

		private double Case1()
		{
			// state_i == 0
			m_z0 = (m_state[m_i + 623] & 0x80000000) | (m_state[m_i + 622] & 0x7FFFFFFF);
			m_z1 = MAT0NEG(-25, m_state[m_i]) ^ MAT0POS(27, m_state[m_i + 70]);
			m_z2 = MAT3POS(9, m_state[m_i + 179]) ^ MAT0POS(1, m_state[m_i + 449]);
			m_state[m_i] = m_z1 ^ m_z2;
			m_state[m_i + 623] = m_z0 ^ MAT0NEG(-9, m_z1) ^ MAT0NEG(-21, m_z2) ^ MAT0POS(21, m_state[m_i]);
			m_i = 623;
			m_case = 3;
			if (TEMPERING)
			{
				 m_y = m_state[m_i] ^ ((m_state[m_i] << 7) & 0xE46E1700);
				 m_y = m_y ^ ((m_y << 15) & 0x9B868000);
				return m_y * FACT;
			}
			else
				return m_state[m_i] * FACT;
		
		}

		private double Case2()
		{
			// state_i == 1
			m_z0 = (m_state[m_i - 1] & 0x80000000) | (m_state[m_i + 622] & 0x7FFFFFFF);
			m_z1 = MAT0NEG(-25, m_state[m_i]) ^ MAT0POS(27, m_state[m_i + 70]);
			m_z2 = MAT3POS(9, m_state[m_i + 179]) ^ MAT0POS(1, m_state[m_i + 449]);
			m_state[m_i] = m_z1 ^ m_z2;
			m_state[m_i - 1] = m_z0 ^ MAT0NEG(-9, m_z1) ^ MAT0NEG(-21, m_z2) ^ MAT0POS(21, m_state[m_i]);
			m_i = 0;
			m_case = 1;
			if (TEMPERING)
			{
				 m_y = m_state[m_i] ^ ((m_state[m_i] << 7) & 0xE46E1700);
				 m_y = m_y ^ ((m_y << 15) & 0x9B868000);
				return m_y * FACT;
			}
			else
				return m_state[m_i] * FACT;
		
		}

		private double Case3()
		{
			// state_i+M1 >= R
			m_z0 = (m_state[m_i - 1] & 0x80000000) | (m_state[m_i - 2] & 0x7FFFFFFF);
			m_z1 = MAT0NEG(-25, m_state[m_i]) ^ MAT0POS(27, m_state[m_i - 554]);
			m_z2 = MAT3POS(9, m_state[m_i - 445]) ^ MAT0POS(1, m_state[m_i - 175]);
			m_state[m_i] = m_z1 ^ m_z2;
			m_state[m_i - 1] = m_z0 ^ MAT0NEG(-9, m_z1) ^ MAT0NEG(-21, m_z2) ^ MAT0POS(21, m_state[m_i]);
			m_i--;
			if (m_i < 554) m_case = 5;
			if (TEMPERING)
			{
				 m_y = m_state[m_i] ^ ((m_state[m_i] << 7) & 0xE46E1700);
				 m_y = m_y ^ ((m_y << 15) & 0x9B868000);
				return m_y * FACT;
			}
			else
				return m_state[m_i] * FACT;
		
		}

		private double Case4()
		{
			// state_i+M3 >= R
			m_z0 = (m_state[m_i - 1] & 0x80000000) | (m_state[m_i - 2] & 0x7FFFFFFF);
			m_z1 = MAT0NEG(-25, m_state[m_i]) ^ MAT0POS(27, m_state[m_i + 70]);
			m_z2 = MAT3POS(9, m_state[m_i + 179]) ^ MAT0POS(1, m_state[m_i - 175]);
			m_state[m_i] = m_z1 ^ m_z2;
			m_state[m_i - 1] = m_z0 ^ MAT0NEG(-9, m_z1) ^ MAT0NEG(-21, m_z2) ^ MAT0POS(21, m_state[m_i]);
			m_i--;
			if (m_i < 175) m_case = 6;

			if (TEMPERING)
			{
				 m_y = m_state[m_i] ^ ((m_state[m_i] << 7) & 0xE46E1700);
				 m_y = m_y ^ ((m_y << 15) & 0x9B868000);
				return m_y * FACT;
			}
			else
				return m_state[m_i] * FACT;
		
		}

		private double Case5()
		{
			// state_i+M2 >= R
			m_z0 = (m_state[m_i - 1] & 0x80000000) | (m_state[m_i - 2] & 0x7FFFFFFF);
			m_z1 = MAT0NEG(-25, m_state[m_i]) ^ MAT0POS(27, m_state[m_i + 70]);
			m_z2 = MAT3POS(9, m_state[m_i - 445]) ^ MAT0POS(1, m_state[m_i - 175]);
			m_state[m_i] = m_z1 ^ m_z2;
			m_state[m_i - 1] = m_z0 ^ MAT0NEG(-9, m_z1) ^ MAT0NEG(-21, m_z2) ^ MAT0POS(21, m_state[m_i]);
			m_i--;
			if (m_i < 445) m_case = 4;

			if (TEMPERING)
			{
				 m_y = m_state[m_i] ^ ((m_state[m_i] << 7) & 0xE46E1700);
				 m_y = m_y ^ ((m_y << 15) & 0x9B868000);
				return m_y * FACT;
			}
			else
				return m_state[m_i] * FACT;
		}

		private double Case6()
		{
			// 2 <= state_i <= (R - M3 - 1)
			m_z0 = (m_state[m_i - 1] & 0x80000000) | (m_state[m_i - 2] & 0x7FFFFFFF);
			m_z1 = MAT0NEG(-25, m_state[m_i]) ^ MAT0POS(27, m_state[m_i + 70]);
			m_z2 = MAT3POS(9, m_state[m_i + 179]) ^ MAT0POS(1, m_state[m_i + 449]);
			m_state[m_i] = m_z1 ^ m_z2;
			m_state[m_i - 1] = m_z0 ^ MAT0NEG(-9, m_z1) ^ MAT0NEG(-21, m_z2) ^ MAT0POS(21, m_state[m_i]);
			m_i--;
			if (m_i == 1) m_case = 2;

			if (TEMPERING)
			{
				 m_y = m_state[m_i] ^ ((m_state[m_i] << 7) & 0xE46E1700);
				 m_y = m_y ^ ((m_y << 15) & 0x9B868000);
				return m_y * FACT;
			}
			else
				return m_state[m_i] * FACT;
		}
	}
}

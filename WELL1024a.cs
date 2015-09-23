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
	public class WELL1024a : Random
	{
		private static double FACT = 2.32830643653869628906e-10;
		private uint	m_i;
		private uint[]	m_state = new uint[32];
		private uint	m_z0;
		private uint	m_z1;
		private uint	m_z2;
		
		public WELL1024a() : this(DateTime.Now.Millisecond)
		{ }
		public WELL1024a(int seed)
		{
			var rnd = new Random(seed);
			for (int i = 0; i < 32; ++i)
				 m_state[i] = unchecked((uint)rnd.Next());
		}
		public WELL1024a(uint[] seed)
		{
			if (seed == null)
				throw new ArgumentNullException("The seed must not be null.");
			if (seed.Length != 32)
				throw new ArgumentException("Length of seed must be 32.");

			for (int i = 0; i < 32; ++i)
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

		protected override double Sample()
		{
			m_z0 = m_state[(m_i + 31) & 0x1F];
			m_z1 = m_state[m_i] ^ MAT0POS(8, m_state[(m_i + 13) & 0x1F]);
			m_z2 = MAT0NEG(-19, m_state[(m_i + 9) & 0x1F]) ^ MAT0NEG(-14, m_state[(m_i + 5) & 0x1F]);
			m_state[m_i] =  m_z1 ^  m_z2;
			m_state[(m_i + 31) & 0x1F] = MAT0NEG(-11,  m_z0) ^ MAT0NEG(-7,  m_z1) ^ MAT0NEG(-13,  m_z2);
			m_i = (m_i + 31) & 0x1F;
			return ((double)m_state[m_i] * FACT);
		}
	}
}

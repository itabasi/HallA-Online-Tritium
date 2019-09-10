

void scale_calc(){





  
  
  TVector3 L_v, R_v, B_v;
  L_v.SetMagThetaPhi( Lp[0], Lth[0], -Lph[0] );
  R_v.SetMagThetaPhi( Rp[0], Rth[0], -Rph[0] );
  B_v.SetMagThetaPhi( sqrt(Ee*Ee-Me*Me), 0, 0 );
  L_v.RotateZ( -13.2 / 180. * PI );
  R_v.RotateZ(  13.2 / 180. * PI );
  mm = sqrt( (Ee + Mp - Ee_ - Ek)*(Ee + Mp - Ee_ - Ek)
                  - (B_v - L_v - R_v)*(B_v - L_v - R_v) );
  




}

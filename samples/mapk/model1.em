Stepper ODEStepper( DE1 ){}
#Stepper DiscreteEventStepper( NR1 ){}


@{MAIN_STEPPER='DE1'}
@{MAIN_PROCESS='MassActionFluxProcess'}
@#{MAIN_STEPPER='NR1'}
@#{MAIN_PROCESS='GillespieProcess'}

@{
VOL = 1e-15
#VOL = 1e-14
N_A = 6.02214e+23
}

@{
def C2N( conc ):
    num = N_A * VOL * conc
    print round( num )
}


System System( / )
{
        StepperID       @MAIN_STEPPER;

        Variable Variable( SIZE )
        {
                Value   @(VOL);
        }

        Variable Variable( K )
        {
                Value   @{C2N( 200e-9 )};
                # SS: 67e-9
        }

        Variable Variable( KK )
        {
                Value   @{C2N( 50e-9 )};
                # SS: 31e-9
        }
        
        
        Variable Variable( P )
        {
                Value   @{C2N( 50e-9 )};
                # SS: 31e-8
        }

        Variable Variable( Kp )
        {
                Value   0;
                # SS: 27e-9
        }

        Variable Variable( Kpp )
        {
                Value   0;
                # SS: 67e-9
        }

        Variable Variable( K_KK )
        {
                Value   0;
                # 16e-9
        }

        Variable Variable( Kp_KK )
        {
                Value   0;
                # 16e-9
        }

        Variable Variable( Kpp_KK )
        {
                Value   0;      
                # 0
                
        }

        Variable Variable( Kpp_P )
        {
                Value   0;
                # 17e-9       
        }

        Variable Variable( Kp_P )
        {
                Value   0;
                # 1.7e-9       
        }

        Process @(MAIN_PROCESS)( R1 )
        {
                VariableReferenceList   [ _ :.:K      -1 ] 
                                        [ _ :.:KK  -1 ]
                                        [ _ :.:K_KK  1];
                k       @( 0.02 * 1e9 );
        }

        Process @(MAIN_PROCESS)( R2 )
        {
                VariableReferenceList   [ _ :.:K_KK -1 ]
                                        [ _ :.:K       1 ] 
                                        [ _ :.:KK   1 ];
                k       1;
        }

        Process @(MAIN_PROCESS)( R3 )
        {
                VariableReferenceList   [ _ :.:K_KK -1 ]
                                        [ _ :.:Kp      1 ] 
                                        [ _ :.:KK   1 ];
                k       1.5;
        }


        Process @(MAIN_PROCESS)( R4 )
        {
                VariableReferenceList   [ _ :.:Kp       -1 ]
                                        [ _ :.:KK     -1 ] 
                                        [ _ :.:Kp_KK   1 ];
                k       @( 0.032 * 1e9 );
        }

        Process @(MAIN_PROCESS)( R5 )
        {
                VariableReferenceList   [ _ :.:Kp_KK  -1 ]
                                        [ _ :.:Kp         1 ]
                                        [ _ :.:KK      1 ];
                k       1;
        }

        Process @(MAIN_PROCESS)( R6 )
        {
                VariableReferenceList   [ _ :.:Kp_KK  -1 ]
                                        [ _ :.:Kpp        1 ]
                                        [ _ :.:KK      1 ]; 
                k       15;
        }

        Process @(MAIN_PROCESS)( R7 )
        {
                VariableReferenceList   [ _ :.:Kpp       -1 ]
                                        [ _ :.:P       -1 ]
                                        [ _ :.:Kpp_P    1 ];
                k       @( 0.02 * 1e9 );
        }

        Process @(MAIN_PROCESS)( R8 )
        {
                VariableReferenceList   [ _ :.:Kpp_P   -1 ] 
                                        [ _ :.:Kpp        1 ]
                                        [ _ :.:P        1 ];
                k       1;
        }

        Process @(MAIN_PROCESS)( R9 )
        {
                VariableReferenceList   [ _ :.:Kpp_P  -1 ]
                                        [ _ :.:Kp        1 ]
                                        [ _ :.:P      1 ]; 
                k       1.5;
        }

        Process @(MAIN_PROCESS)( R10 )
        {
                VariableReferenceList   [ _ :.:Kp       -1 ]
                                        [ _ :.:P       -1 ]
                                        [ _ :.:Kp_P    1 ];
                k       @( 0.032 * 1e9 );
        }

        Process @(MAIN_PROCESS)( R11 )
        {
                VariableReferenceList   [ _ :.:Kp_P       -1 ]
                                        [ _ :.:Kp       1 ]
                                        [ _ :.:P    1 ];
                k      1;
        }

        Process @(MAIN_PROCESS)( R12 )
        {
                VariableReferenceList   [ _ :.:Kp_P   -1 ] 
                                        [ _ :.:K        1 ]
                                        [ _ :.:P        1 ];
                k       15;
        }

        

        

}


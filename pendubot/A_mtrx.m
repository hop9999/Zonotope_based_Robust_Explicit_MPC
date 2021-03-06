function A_mtrx = A_mtrx(in1,in2,u1,in4,me1,le1,me2,le2)
%A_MTRX
%    A_MTRX = A_MTRX(IN1,IN2,U1,IN4,ME1,LE1,ME2,LE2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    21-Feb-2021 16:45:44

dq1 = in2(1,:);
dq2 = in2(2,:);
p1 = in4(1,:);
p2 = in4(2,:);
p3 = in4(3,:);
p4 = in4(4,:);
p5 = in4(5,:);
p6 = in4(6,:);
p7 = in4(7,:);
p8 = in4(8,:);
p9 = in4(9,:);
q1 = in1(1,:);
q2 = in1(2,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = tanh(dq1);
t5 = tanh(dq2);
t6 = sin(q1);
t7 = sin(q2);
t8 = le2.*me2;
t9 = me2.*p3;
t10 = q1+q2;
t11 = dq1.^2;
t12 = dq2.^2;
t13 = le1.^2;
t14 = le2.^2;
t15 = le2.^3;
t16 = me2.^2;
t17 = p2.*4.0;
t18 = p2.^2;
t19 = p3.^2;
t20 = q2.*2.0;
t29 = -q2;
t31 = p1.*p3.*1.6e+1;
t32 = p1.*p3.*3.2e+1;
t21 = cos(t20);
t22 = t3.^2;
t23 = t4.^2;
t24 = t5.^2;
t25 = t9.*2.0;
t26 = sin(t20);
t27 = cos(t10);
t28 = sin(t10);
t33 = q2+t10;
t35 = t18.*1.6e+1;
t36 = t19.*1.6e+1;
t37 = t19.*3.2e+1;
t38 = -t32;
t39 = q1+t29;
t40 = p2.*t8.*8.0;
t43 = t8+t17;
t46 = t8.^2;
t47 = me1.*p3.*t13.*1.6e+1;
t48 = le2.*p1.*t8.*1.6e+1;
t49 = le2.*p3.*t8.*1.6e+1;
t50 = me1.*p3.*t13.*3.2e+1;
t51 = le2.*p1.*t8.*3.2e+1;
t52 = le2.*p3.*t8.*3.2e+1;
t53 = t6.*t9.*1.96e+2;
t57 = p3.*p4.*t6.*7.84e+2;
t62 = dq1.*p2.*p3.*t7.*4.0e+3;
t63 = dq2.*p2.*p3.*t7.*4.0e+3;
t64 = le1.*me1.*p3.*t6.*7.84e+2;
t70 = dq2.*p2.*p3.*t7.*1.6e+4;
t71 = le2.*me1.*t8.*t13.*1.6e+1;
t72 = le2.*me1.*t8.*t13.*3.2e+1;
t74 = dq1.*p3.*t7.*t8.*1.0e+3;
t75 = dq2.*p3.*t7.*t8.*1.0e+3;
t76 = dq2.*p3.*t7.*t8.*4.0e+3;
t79 = le2.*p4.*t6.*t8.*7.84e+2;
t82 = dq1.*le2.*p2.*t7.*t8.*4.0e+3;
t83 = dq2.*le2.*p2.*t7.*t8.*4.0e+3;
t84 = le1.*le2.*me1.*t6.*t8.*7.84e+2;
t85 = dq2.*p2.*t3.*t7.*t8.*8.0e+3;
t93 = dq2.*t3.*t7.*t18.*1.6e+4;
t94 = dq1.*le2.*p2.*t7.*t8.*1.6e+4;
t95 = dq2.*le2.*p2.*t7.*t8.*1.6e+4;
t30 = -t25;
t34 = sin(t33);
t41 = sin(t39);
t42 = -t36;
t44 = t23-1.0;
t45 = t24-1.0;
t54 = -t49;
t55 = -t50;
t56 = -t51;
t58 = t43.^2;
t59 = -t46;
t60 = t21.*t35;
t61 = t22.*t35;
t65 = t21.*t40;
t66 = t22.*t40;
t67 = t21.*t46;
t68 = t22.*t46;
t69 = t18.*t22.*-1.6e+1;
t73 = p2.*t8.*t22.*-8.0;
t77 = -t72;
t78 = t6.*t46.*1.96e+2;
t81 = p5.*t3.*t8.*t28.*1.96e+2;
t86 = p2.*p5.*t3.*t28.*7.84e+2;
t88 = p2.*t3.*t8.*t28.*7.84e+2;
t89 = dq1.*le2.*t7.*t46.*1.0e+3;
t90 = dq2.*le2.*t7.*t46.*1.0e+3;
t91 = dq1.*le2.*t7.*t46.*4.0e+3;
t92 = dq2.*le2.*t7.*t46.*4.0e+3;
t98 = dq2.*t3.*t7.*t46.*1.0e+3;
t99 = t3.*t28.*t46.*1.96e+2;
t80 = t22.*t59;
t87 = -t81;
t96 = -t86;
t97 = -t88;
t100 = -t99;
t104 = t30+t35+t37+t38+t40+t52+t55+t56+t59+t60+t65+t67+t77;
t101 = t9+t31+t42+t46+t47+t48+t54+t69+t71+t73+t80;
t105 = 1.0./t104;
t102 = 1.0./t101;
t103 = t102.^2;
A_mtrx = reshape([0.0,0.0,(t102.*(t53+t57+t64+t78+t79+t84+t87+t96+t97+t100))./5.0,t102.*(t53+t57+t64+t78+t79+t84+t87+t96+t97+t100-me2.*p5.*t28.*4.9e+1-p1.*p5.*t28.*7.84e+2+p3.*p5.*t28.*7.84e+2-me2.*t8.*t28.*4.9e+1-p1.*t8.*t28.*7.84e+2+p3.*t8.*t28.*7.84e+2+me2.*p2.*t3.*t6.*1.96e+2-me1.*p5.*t13.*t28.*7.84e+2+p2.*p4.*t3.*t6.*7.84e+2+me2.*t3.*t6.*t8.*4.9e+1-me1.*t8.*t13.*t28.*7.84e+2+p4.*t3.*t6.*t8.*1.96e+2+le1.*me1.*p2.*t3.*t6.*7.84e+2+le1.*me1.*t3.*t6.*t8.*1.96e+2).*(-1.0./5.0),0.0,0.0,-t105.*(p2.*(8.0./5.0)+t8.*(2.0./5.0)).*(p5.*t34.*-1.96e+2-t8.*t34.*1.96e+2-dq2.*p7.*t7.*2.0e+1+p3.*t3.*t11.*2.0e+1+p3.*t3.*t12.*2.0e+1-p9.*t5.*t7.*2.0e+1+p2.*t11.*t21.*2.0e+1+t8.*t11.*t21.*5.0+dq1.*dq2.*p3.*t3.*4.0e+1+le2.*t3.*t8.*t11.*2.0e+1+le2.*t3.*t8.*t12.*2.0e+1+dq1.*dq2.*le2.*t3.*t8.*4.0e+1)-(t26.*t58.*t103.*(dq2.*t62+dq2.*t74+dq2.*t82+dq2.*t89+p3.*u1.*2.26e+2-t2.*t9.*4.9e+3-t2.*t46.*4.9e+3-dq1.*p3.*p6.*2.0e+3+dq2.*p3.*p7.*2.0e+3-p3.*p4.*t2.*1.96e+4-p3.*p8.*t4.*2.0e+3+p3.*p9.*t5.*2.0e+3+le2.*t8.*u1.*2.26e+2+t11.*t18.*t26.*1.0e+3+t3.*t27.*t46.*4.9e+3+t11.*t26.*t46.*(1.25e+2./2.0)-dq1.*le2.*p6.*t8.*2.0e+3+dq2.*le2.*p7.*t8.*2.0e+3+dq2.*p2.*p7.*t3.*2.0e+3+dq2.*p7.*t3.*t8.*5.0e+2-le1.*me1.*p3.*t2.*1.96e+4-le2.*p4.*t2.*t8.*1.96e+4-le2.*p8.*t4.*t8.*2.0e+3+le2.*p9.*t5.*t8.*2.0e+3+p2.*p9.*t3.*t5.*2.0e+3+p2.*p3.*t7.*t11.*2.0e+3+p2.*p3.*t7.*t12.*2.0e+3+p2.*p5.*t3.*t27.*1.96e+4+le2.*t7.*t11.*t46.*5.0e+2+le2.*t7.*t12.*t46.*5.0e+2+p9.*t3.*t5.*t8.*5.0e+2+p3.*t7.*t8.*t11.*5.0e+2+p3.*t7.*t8.*t12.*5.0e+2+p2.*t3.*t8.*t27.*1.96e+4+p5.*t3.*t8.*t27.*4.9e+3+p2.*t8.*t11.*t26.*5.0e+2-le1.*le2.*me1.*t2.*t8.*1.96e+4+le2.*p2.*t7.*t8.*t11.*2.0e+3+le2.*p2.*t7.*t8.*t12.*2.0e+3))./1.25e+2,(t105.*(t11.*t67.*1.0e+3+t12.*t67.*5.0e+2-t34.*t46.*1.96e+4+dq1.*dq2.*t67.*1.0e+3+me2.*p2.*t28.*9.8e+3-me2.*p5.*t28.*4.9e+3-me2.*p2.*t41.*9.8e+3-p1.*p5.*t28.*7.84e+4+p2.*p4.*t28.*3.92e+4+p3.*p5.*t28.*7.84e+4-p2.*p5.*t34.*7.84e+4-p2.*p4.*t41.*3.92e+4-me2.*t8.*t28.*2.45e+3-me2.*t8.*t41.*2.45e+3-p1.*t8.*t28.*7.84e+4+p3.*t8.*t28.*7.84e+4+p4.*t8.*t28.*9.8e+3-p2.*t8.*t34.*7.84e+4-p5.*t8.*t34.*1.96e+4-p4.*t8.*t41.*9.8e+3-p2.*t7.*u1.*9.04e+2+t11.*t18.*t21.*1.6e+4+t12.*t18.*t21.*8.0e+3-t7.*t8.*u1.*2.26e+2+dq1.*dq2.*t18.*t21.*1.6e+4+dq1.*p2.*p6.*t7.*8.0e+3-dq2.*p2.*p7.*t7.*1.6e+4+dq1.*p6.*t7.*t8.*2.0e+3-dq2.*p7.*t7.*t8.*4.0e+3+le1.*me1.*p2.*t28.*3.92e+4-le1.*me1.*p2.*t41.*3.92e+4+le1.*me1.*t8.*t28.*9.8e+3-le1.*me1.*t8.*t41.*9.8e+3+me2.*p2.*t3.*t11.*5.0e+2-me1.*p5.*t13.*t28.*7.84e+4+p1.*p2.*t3.*t11.*8.0e+3+p2.*p3.*t3.*t12.*8.0e+3+p2.*p8.*t4.*t7.*8.0e+3-p2.*p9.*t5.*t7.*1.6e+4+le2.*t3.*t11.*t46.*2.0e+3+le2.*t3.*t12.*t46.*2.0e+3+me2.*t3.*t8.*t11.*1.25e+2-me1.*t8.*t13.*t28.*7.84e+4+p1.*t3.*t8.*t11.*2.0e+3+p3.*t3.*t8.*t12.*2.0e+3+p8.*t4.*t7.*t8.*2.0e+3-p9.*t5.*t7.*t8.*4.0e+3+p2.*t8.*t11.*t21.*8.0e+3+p2.*t8.*t12.*t21.*4.0e+3+dq1.*dq2.*p2.*p3.*t3.*1.6e+4+dq1.*dq2.*le2.*t3.*t46.*4.0e+3+dq1.*dq2.*p3.*t3.*t8.*4.0e+3+dq1.*dq2.*p2.*t8.*t21.*8.0e+3+le2.*p2.*t3.*t8.*t11.*8.0e+3+le2.*p2.*t3.*t8.*t12.*8.0e+3+me1.*p2.*t3.*t11.*t13.*8.0e+3+me1.*t3.*t8.*t11.*t13.*2.0e+3+dq1.*dq2.*le2.*p2.*t3.*t8.*1.6e+4))./2.5e+2+(t26.*t58.*t103.*(dq1.*t70+dq1.*t76+dq2.*t91+dq2.*t94+p3.*u1.*9.04e+2-t2.*t9.*1.96e+4-t2.*t46.*1.96e+4+dq2.*me2.*p7.*5.0e+2-dq1.*p3.*p6.*8.0e+3+dq2.*p1.*p7.*8.0e+3+me2.*p9.*t5.*5.0e+2+me2.*p5.*t27.*4.9e+3-p3.*p4.*t2.*7.84e+4+p1.*p9.*t5.*8.0e+3-p3.*p8.*t4.*8.0e+3+p1.*p5.*t27.*7.84e+4-p3.*p5.*t27.*7.84e+4+le2.*t8.*u1.*9.04e+2+me2.*t8.*t27.*4.9e+3+p1.*t8.*t27.*7.84e+4-p3.*t8.*t27.*7.84e+4+p2.*t3.*u1.*9.04e+2+t11.*t18.*t26.*8.0e+3+t12.*t18.*t26.*4.0e+3+t3.*t27.*t46.*1.96e+4+t11.*t26.*t46.*5.0e+2+t12.*t26.*t46.*2.5e+2+t3.*t8.*u1.*2.26e+2+dq1.*dq2.*t18.*t26.*8.0e+3+dq1.*dq2.*t26.*t46.*5.0e+2-dq1.*le2.*p6.*t8.*8.0e+3+dq2.*le2.*p7.*t8.*8.0e+3+dq2.*me1.*p7.*t13.*8.0e+3-dq1.*p2.*p6.*t3.*8.0e+3+dq2.*p2.*p7.*t3.*1.6e+4-dq1.*p6.*t3.*t8.*2.0e+3+dq2.*p7.*t3.*t8.*4.0e+3-le1.*me1.*p3.*t2.*7.84e+4-le2.*p4.*t2.*t8.*7.84e+4-le2.*p8.*t4.*t8.*8.0e+3+le2.*p9.*t5.*t8.*8.0e+3-me2.*p2.*t2.*t3.*1.96e+4+me2.*p2.*t7.*t11.*5.0e+2+me1.*p9.*t5.*t13.*8.0e+3+me1.*p5.*t13.*t27.*7.84e+4-p2.*p4.*t2.*t3.*7.84e+4-p2.*p8.*t3.*t4.*8.0e+3+p2.*p9.*t3.*t5.*1.6e+4+p1.*p2.*t7.*t11.*8.0e+3+p2.*p3.*t7.*t12.*8.0e+3+p2.*p5.*t3.*t27.*7.84e+4+le2.*t7.*t11.*t46.*2.0e+3+le2.*t7.*t12.*t46.*2.0e+3-me2.*t2.*t3.*t8.*4.9e+3+me2.*t7.*t8.*t11.*1.25e+2+me1.*t8.*t13.*t27.*7.84e+4-p4.*t2.*t3.*t8.*1.96e+4-p8.*t3.*t4.*t8.*2.0e+3+p9.*t3.*t5.*t8.*4.0e+3+p1.*t7.*t8.*t11.*2.0e+3+p3.*t7.*t8.*t12.*2.0e+3+p2.*t3.*t8.*t27.*7.84e+4+p5.*t3.*t8.*t27.*1.96e+4+p2.*t8.*t11.*t26.*4.0e+3+p2.*t8.*t12.*t26.*2.0e+3+dq1.*dq2.*p2.*t8.*t26.*4.0e+3-le1.*le2.*me1.*t2.*t8.*7.84e+4-le1.*me1.*p2.*t2.*t3.*7.84e+4-le1.*me1.*t2.*t3.*t8.*1.96e+4+le2.*p2.*t7.*t8.*t11.*8.0e+3+le2.*p2.*t7.*t8.*t12.*8.0e+3+me1.*p2.*t7.*t11.*t13.*8.0e+3+me1.*t7.*t8.*t11.*t13.*2.0e+3))./5.0e+2,1.0,0.0,(t102.*(t62+t63+t74+t75+t82+t83+t89+t90-p3.*p6.*2.0e+3-le2.*p6.*t8.*2.0e+3+p3.*p8.*t44.*2.0e+3+dq1.*t3.*t7.*t18.*4.0e+3+dq1.*t3.*t7.*t46.*2.5e+2+le2.*p8.*t8.*t44.*2.0e+3+dq1.*p2.*t3.*t7.*t8.*2.0e+3))./1.25e+2,t102.*(t70+t76+t85+t91+t92+t93+t94+t95+t98-p3.*p6.*8.0e+3-le2.*p6.*t8.*8.0e+3-p2.*p6.*t3.*8.0e+3+p3.*p8.*t44.*8.0e+3-p6.*t3.*t8.*2.0e+3+dq1.*me2.*p2.*t7.*1.0e+3+dq1.*p1.*p2.*t7.*1.6e+4+dq1.*me2.*t7.*t8.*2.5e+2+dq1.*p1.*t7.*t8.*4.0e+3+dq1.*t3.*t7.*t18.*3.2e+4+dq1.*t3.*t7.*t46.*2.0e+3+le2.*p8.*t8.*t44.*8.0e+3+p2.*p8.*t3.*t44.*8.0e+3+p8.*t3.*t8.*t44.*2.0e+3+dq1.*me1.*p2.*t7.*t13.*1.6e+4+dq1.*me1.*t7.*t8.*t13.*4.0e+3+dq1.*p2.*t3.*t7.*t8.*1.6e+4).*(-1.0./5.0e+2),0.0,1.0,(t102.*(t62+t63+t74+t75+t82+t83+t89+t90+p3.*p7.*2.0e+3+le2.*p7.*t8.*2.0e+3+p2.*p7.*t3.*2.0e+3-p3.*p9.*t45.*2.0e+3+p7.*t3.*t8.*5.0e+2-le2.*p9.*t8.*t45.*2.0e+3-p2.*p9.*t3.*t45.*2.0e+3-p9.*t3.*t8.*t45.*5.0e+2))./1.25e+2,t102.*(t70+t76+t85+t91+t92+t93+t94+t95+t98+me2.*p7.*5.0e+2+p1.*p7.*8.0e+3+le2.*p7.*t8.*8.0e+3+me1.*p7.*t13.*8.0e+3-me2.*p9.*t45.*5.0e+2+p2.*p7.*t3.*1.6e+4-p1.*p9.*t45.*8.0e+3+p7.*t3.*t8.*4.0e+3+dq1.*p2.*p3.*t7.*1.6e+4+dq1.*p3.*t7.*t8.*4.0e+3+dq1.*t3.*t7.*t18.*1.6e+4+dq1.*t3.*t7.*t46.*1.0e+3-le2.*p9.*t8.*t45.*8.0e+3-me1.*p9.*t13.*t45.*8.0e+3-p2.*p9.*t3.*t45.*1.6e+4-p9.*t3.*t8.*t45.*4.0e+3+dq1.*p2.*t3.*t7.*t8.*8.0e+3).*(-1.0./5.0e+2)],[4,4]);

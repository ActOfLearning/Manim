from manimlib.imports import *
import numpy as np

def Circle_Inversion(obj, about_point = ORIGIN, inversion_radius = 1.0):
    def inverse_of_a_point(pt0):
        norm_of_pt = np.dot(pt0, pt0)
        if norm_of_pt == 0:
            norm_of_pt += 0.0001
        return pt0 / norm_of_pt
    def inverse_of_a_circle(obj):
        center_of_the_circle = obj.get_center()
        radius_of_the_circle = np.linalg.norm(center_of_the_circle - obj.get_start())
        inversion_factor = np.dot(center_of_the_circle, center_of_the_circle) - radius_of_the_circle ** 2
        if inversion_factor == 0:
            res_object = Line(10 * UP, 10 * DOWN)
            res_object.move_to(RIGHT * (1 / 2 / radius_of_the_circle))
            rot_angle = np.angle(center_of_the_circle[0] + center_of_the_circle[1] * 1j)
            res_object.rotate_about_origin(rot_angle)
            return res_object
        inversion_factor = 1.0 / inversion_factor
        return Circle(radius = abs(inversion_factor) * radius_of_the_circle).move_to(inversion_factor * center_of_the_circle).set_color(TEAL)
    def circum_center(pt0):
        pt = pt0[:]
        side_a, side_b, side_c = pt[2] - pt[1], pt[0] - pt[2], pt[1] - pt[0]
        side_a, side_b, side_c = np.dot(side_a, side_a), np.dot(side_b, side_b), np.dot(side_c, side_c)
        circum_x, circum_y, circum_z = side_a * (side_b + side_c - side_a), side_b * (side_c + side_a - side_b), side_c * (side_a + side_b - side_c)
        return (circum_x * pt[0] + circum_y * pt[1] + circum_z * pt[2]) / (circum_x + circum_y + circum_z)
    def inverted_arc(c0, r0, ang0):
        return Arc(ang0[0], ang0[1] - ang0[0] + (2 * PI if ang0[1] < ang0[0] else 0), radius = r0, arc_center = c0)
    def arg_of_a_point(pt0):
        return np.angle(pt0[0] + 1j * pt0[1])
    def inverse_of_a_line(obj):
        line_st, line_en = obj.get_start(), obj.get_end()
        inv_line_start = inverse_of_a_point(line_st)
        inv_line_end = inverse_of_a_point(line_en)
        if abs(np.linalg.det([line_st, line_en, OUT])) >= 0.001:
            circum = circum_center([ORIGIN, inv_line_start, inv_line_end])
            circum_norm = np.linalg.norm(circum)
            temp_mat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            temp_mat = np.dot(np.dot(temp_mat, np.transpose(line_st)), line_en)
            return inverted_arc(circum, circum_norm, [arg_of_a_point(inv_line_start - circum), arg_of_a_point(inv_line_end - circum)] if temp_mat > 0 else [arg_of_a_point(inv_line_end - circum), arg_of_a_point(inv_line_start - circum)])
        else:
            if np.dot(line_st, line_en) > 0:
                return Line(inv_line_start, inv_line_end)
            else:
                return VGroup(Line(inv_line_start, 6 * inv_line_start / (0.1 + np.linalg.norm(inv_line_start))), Line(inv_line_end, 6 * inv_line_end / (0.1 + np.linalg.norm(inv_line_end))))
    def inverse_of_a_polygon(obj):
        polygon_vertices = obj.get_points()
        inverted_points = VGroup()
        for i in range(1, len(polygon_vertices)):
            inverted_points.add(inverse_of_a_line(Line(polygon_vertices[i], polygon_vertices[i - 1])))
        return inverted_points.set_color(TEAL)
    shifted_obj = obj.copy()
    shifted_obj.shift(-about_point)
    shifted_obj.scale(1 / inversion_radius, about_point = ORIGIN)
    if type(obj).__name__ == "Circle":
        res_obj = inverse_of_a_circle(shifted_obj)
    elif type(obj).__name__ == "Dot":
        res_obj = Dot(radius = DEFAULT_DOT_RADIUS / inversion_radius).move_to(inverse_of_a_point(shifted_obj.get_center()))
    elif type(obj).__name__ == "VGroup":
        res_obj = VGroup()
        for objl in shifted_obj:
            res_obj.add(Circle_Inversion(objl))
    else:
        res_obj = inverse_of_a_polygon(shifted_obj)
    res_obj.scale(inversion_radius, about_point = ORIGIN)
    res_obj.shift(about_point)
    return res_obj.set_color(TEAL)

def ParallalLine(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            return Line(p_cent, p_cent + normalized_line)
    return res

def PerpendicularFoot(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            line_start += np.dot(p_cent - line_start, normalized_line) * normalized_line
            return Dot().move_to(line_start)
    return res

def PerpendicularLine(pt0, line0):
    line_start, line_end = line0.get_start_and_end()
    line_unit = (line_end - line_start) / np.linalg.norm(line_end - line_start)
    perp_unit = np.cross(line_unit, OUT)
    perp_foot = PerpendicularFoot(pt0, line0).get_center()
    return Line(perp_foot, perp_unit + perp_foot)

def DottoCircle(obj1):
    cen1, rad1 = ORIGIN, 0.0
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1 = obj1.get_center()
        if type(obj1).__name__ == "Dot":
            rad1 = 0.0
        else:
            rad1 = np.linalg.norm(cen1 - obj1.get_start())
    return cen1, rad1

def RadicalAxis(obj1, obj2):
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1, rad1 = DottoCircle(obj1)
    if type(obj2).__name__ == "Dot" or type(obj2).__name__ == "Circle":
        cen2, rad2 = DottoCircle(obj2)
    dist_cent = np.linalg.norm(cen1 - cen2)
    if dist_cent > 0:
        res = rad1 ** 2 - rad2 ** 2
        res = dist_cent + res / dist_cent
        res /= 2
        res /= dist_cent
        res = cen1 + res * (cen2 - cen1)
        return PerpendicularLine(Dot().move_to(res), Line(cen1, cen2))

def LineSegIntersection(obj1, obj2):
    res_obj = VGroup()
    obj1_start, obj1_end = obj1.get_start_and_end()
    a = obj1_end - obj1_start
    obj2_start, obj2_end = obj2.get_start_and_end()
    b = obj2_end - obj2_start
    c = obj2_start - obj1_start
    acrossb = np.cross(a, b)
    bcrossa = np.cross(b, a)
    acrossc = np.cross(a, c)
    ccrossb = np.cross(c, b)
    acrossb_dot = np.dot(acrossb, acrossb)
    bcrossa_dot = np.dot(bcrossa, bcrossa)
    if acrossb_dot > 0:
        s = np.dot(ccrossb, acrossb) / acrossb_dot
        t = np.dot(acrossc, bcrossa) / bcrossa_dot
        if 0 <= s <= 1 and 0 <= t <= 1:
            res_obj.add(Dot(obj1_start + s * a))
    return res_obj

def RayIntersection(obj1, obj2):
    #needs more work
    res_obj = VGroup()
    obj1_start, obj1_end = obj1.get_start_and_end()
    ab = obj1_end - obj1_start
    ab_unit = ab / np.linalg.norm(ab)
    obj2_start, obj2_end = obj2.get_start_and_end()
    cd = obj2_end - obj2_start
    cd_unit = cd / np.linalg.norm(cd)
    ac = obj2_start - obj1_start
    v = np.cross(ab_unit, cd_unit)
    if np.linalg.norm(v) > 0:
        v_unit = v / np.linalg.norm(v)
        v_sq = np.dot(v, v)
        s1 = np.dot(ac, np.cross(cd_unit, v_unit)) / v_sq
        s2 = np.dot(ac, np.cross(ab_unit, v_unit)) / v_sq
        print(s1, s2, obj2_start, obj2_end)
        if 0 < s1 < 1 and 0 < s2 < 1:
            return (obj1_start + ab_unit * s1 + obj2_start + cd_unit * s2) / 2
    return res_obj

def Intersection(*args):
    res_obj = VGroup()
    obj1, obj2 = args[0], args[1]
    if type(obj1).__name__ == "Line":
        obj1_start, obj1_end = obj1.get_start_and_end()
        ab = obj1_end - obj1_start
        ab_unit = ab / np.linalg.norm(ab)
        if type(obj2).__name__ == "Line":
            obj2_start, obj2_end = obj2.get_start_and_end()
            cd = obj2_end - obj2_start
            ac = obj2_start - obj1_start
            abcd = np.cross(ab, cd)
            abcd_n = np.dot(abcd, abcd)
            if abcd_n != 0:
                abcd = abcd / abcd_n
                abcd = np.dot(np.cross(ac, cd), abcd)
                abcd = obj1_start + ab * abcd
                res_obj.add(Dot(abcd))
        elif type(obj2).__name__ == "Circle":
            cen = Dot().move_to(obj2.get_center())
            p_foot = PerpendicularFoot(cen, obj1)
            p_dist = np.linalg.norm(cen.get_center() - p_foot.get_center())
            rad = np.linalg.norm(obj2.get_center() - obj2.get_start())
            if p_dist == rad:
                res_obj.add(p_foot)
            elif p_dist < rad:
                chord_dist = rad * rad - p_dist * p_dist
                chord_dist = np.sqrt(chord_dist)
                res_obj.add(Dot().move_to(p_foot.get_center() + chord_dist * ab_unit), Dot().move_to(p_foot.get_center() - chord_dist * ab_unit))
    elif type(obj1).__name__ == "Circle":
        cen1, cen2 = obj1.get_center(), obj2.get_center()
        dist_cen = np.linalg.norm(cen1 - cen2)
        rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
        if abs(dist_cen - rad1 - rad2) < 1e-3:
            res_obj.add(Dot().move_to((rad2 * cen1 + rad1 * cen2) / (rad1 + rad2)))
        elif dist_cen < rad1 + rad2:
            p_foot = RadicalAxis(obj1, obj2)
            #p_foot = PerpendicularLine(p_foot, Line(cen1, cen2))
            return Intersection(p_foot, obj1)
    return res_obj

def Midpoint(*args):
    if type(args[0]).__name__ == "Line":
        line_start, line_end = args[0].get_start_and_end()
        mid_pt = (line_start + line_end) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        mid_pt = (args[0].get_center() + args[1].get_center()) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt

def PerpendicularBisector(*args):
    mid_pt = Midpoint(*args)
    if type(args[0]).__name__ == "Line":
        print("It's a line..")
        return PerpendicularLine(mid_pt, args[0])
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        connec_line = Line(args[0].get_center(), args[1].get_center())
        return PerpendicularLine(mid_pt, connec_line)

def AngularBisector(*args):
    if type(args[0]).__name__ == "Line" and type(args[1]).__name__ == "Line":
        intersec = Intersection(*args)
        if len(list(intersec)) > 0:
            return AngularBisector(Dot().move_to(args[0].get_end()), intersec, Dot().move_to(args[1].get_end()))
        else:
            return VGroup()
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot" and type(args[2]).__name__ == "Dot":
        pt_a = args[0].get_center()
        pt_b = args[1].get_center()
        pt_c = args[2].get_center()
        vec_ba = pt_a - pt_b
        vec_bc = pt_c - pt_b
        return Line(pt_b, vec_ba + vec_bc)

def RadicalCenter(*args):
    if type(args[0]).__name__ == "Circle" and type(args[1]).__name__ == "Circle" and type(args[2]).__name__ == "Circle":
        radic12 = RadicalAxis(args[0], args[1])
        radic23 = RadicalAxis(args[1], args[2])
        return Intersection(radic12, radic23)
    return Dot()

def HomotheticCenter(obj1, obj2):
    cen1, cen2 = obj1.get_center(), obj2.get_center()
    rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
    res_obj = VGroup()
    res_obj.add(Dot().move_to((rad1 * cen2 + rad2 * cen1) / (rad1 + rad2)))
    if rad1 != rad2:
        res_obj.add(Dot().move_to((rad1 * cen2 - rad2 * cen1) / (rad1 - rad2)))
    return res_obj

def CircleTangentLines(obj, pt):
    cen = obj.get_center()
    rad = np.linalg.norm(cen - obj.get_start())
    dist = np.linalg.norm(cen - pt.get_center())
    res_obj = VGroup()
    if dist < rad:
        return res_obj.add(pt)
    temp_c = Circle(arc_center = (cen + pt.get_center()) / 2, radius = dist / 2)
    temp_c = Intersection(temp_c, obj)
    for i in temp_c:
        res_obj.add(Line(pt.get_center(), i.get_center()))
    return res_obj

class Angle(VGroup):

    CONFIG = {
        'radius': 0.5,
        'color': RED,
        'opacity': 0.4,
        'stroke_width': 5,
        'below_180': True,
    }

    def __init__(self, A, O, B, **kwargs):

        VMobject.__init__(self, **kwargs)
        OA, OB = A-O, B-O
        if self.below_180:
            theta = np.angle(complex(*OA[:2])/complex(*OB[:2])) # angle of OB to OA
        else:
            theta = TAU + np.angle(complex(*OA[:2])/complex(*OB[:2]))

        if abs(abs(theta) - PI / 2) <= 1e-3:
            self.add(Elbow(angle = np.angle(complex(*OB[:2])), color = self.color, width = self.radius).shift(O))
        else:
            '''self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=self.radius/2,
                     stroke_width=100 * self.radius, color=self.color, arc_center=O).set_stroke(opacity=self.opacity))'''
            #self.add(Sector(inner_radius = 0, outer_radius = self.radius, start_angle = Line(O, B), angle = theta, color = self.color, arc_center = O, fill_opacity = self.opacity))
            self.add(Sector(inner_radius = 0, outer_radius = self.radius, angle = theta, arc_center = O, color = self.color, fill_opacity = self.opacity, start_angle = Line(O, B).get_angle()))
            self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=self.radius,
                     stroke_width=self.stroke_width, color=self.color, arc_center=O))

class LogSpiral(ParametricFunction):
    CONFIG = {
        "radius"        : 1,
        "alpha"         : 30 * DEGREES,
        "start_val"   : -2 * PI,
        "end_val"     : 2 * PI,
        "step_size"     : 0.001,
        #"density"       : 50 * DEFAULT_POINT_DENSITY_1D,
        "color"         : BLUE,
    }
    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        ParametricFunction.__init__(self, self.pos_func, **kwargs)

    def pos_func(self, t):
        T = self.start_val + t * (self.end_val - self.start_val)
        return self.radius * np.array([
            np.exp(T / np.tan(self.alpha)) * np.cos(T),
            np.exp(T / np.tan(self.alpha)) * np.sin(T),
            0
        ])

class Cycloid(ParametricFunction):
    CONFIG = {
        "point_a"       : ORIGIN,
        "radius"        : 1,
        "start_theta"   : 0,
        "end_theta"     : 2 * np.pi,
        "density"       : 5 * DEFAULT_POINT_DENSITY_1D,
        "color"         : BLUE,
        "inverted"      : False,
        "frac"          : 1.0
    }
    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        ParametricFunction.__init__(self, self.pos_func, **kwargs)

    def pos_func(self, t):
        T = self.start_theta + t * (self.end_theta - self.start_theta)
        return self.point_a + np.array([
            self.radius * T - (-1 if self.inverted else 1) * self.frac * self.radius * np.sin(T),
            self.radius - self.radius * self.frac * np.cos(T),
            0
        ])

class ProblemIntroduction(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(4.5 * DOWN, 4.5 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        polaraxes = VGroup()
        num_lines, rad_inc = 10, 1
        start_ang = 180 * DEGREES / num_lines
        start_rad = rad_inc
        while start_ang < 180 * DEGREES:
            polaraxes.add(Line(10 * LEFT, 10 * RIGHT).rotate(start_ang))
            start_ang += 360 * DEGREES / num_lines
        while start_rad <= 10:
            polaraxes.add(Circle(radius = start_rad))
            start_rad += rad_inc
        polaraxes.set_color(GREY)
        polaraxes.fade(0.875)
        self.play(ShowCreation(polaraxes), run_time = 3)
        self.wait()
        tnplane, tpolaraxes = nplane.copy(), polaraxes.copy()
        #spiang = 45 * DEGREES * (1 + 1 / 2 + 1 / 8)
        spiang = 82.5 * DEGREES
        #spiang = ValueTracker(np.arctan(3.25))
        spi = LogSpiral(start_val = -25, end_val = 25, alpha = spiang)
        self.play(ShowCreation(spi), run_time = 7)
        self.wait()
        line, rad_line = Line(ORIGIN, 10 * RIGHT).fade(0.75), Line(ORIGIN, RIGHT).set_color(GREEN)
        tracingdot = Dot(rad_line.get_end(), color = RED)
        path = VMobject(color = RED)
        path.set_points_as_corners([tracingdot.get_center(), tracingdot.get_center() + 0.001 * DOWN])
        self.play(FadeOut(spi))
        self.wait()
        alp = ValueTracker(0)
        omega = 1
        linesgp = VGroup(line, rad_line)
        path_gp = VGroup(tracingdot, path)
        def lineupdater(obj):
            li, rli = obj
            li.become(Line(ORIGIN, 10 * RIGHT))
            li.rotate_about_origin(omega * alp.get_value())
            li.fade(0.75)
            leng = np.exp(omega * alp.get_value() / np.tan(spiang))
            temprli = Line(ORIGIN, leng * RIGHT).set_color(GREEN)
            temprli.rotate_about_origin(omega * alp.get_value())
            rli.become(temprli)
        def pathupdater(obj):
            ddot, ppath = obj
            ddot.move_to(rad_line.get_end())
            oldpath = ppath.copy()
            oldpath.append_vectorized_mobject(Line(oldpath.points[-1], ddot.get_center()))
            oldpath.make_smooth()
            ppath.become(oldpath)
        self.play(Write(linesgp), Write(path_gp))
        linesgp.add_updater(lineupdater)
        path_gp.add_updater(pathupdater)
        self.add(linesgp, path_gp)
        velocityrelation = TextMobject("The velocity of the tracing point along the line is directly proportional to its from the origin (or the pole) while the line itself is rotating at an uniform angular velocity.").scale(3 / 4)
        velbrec = BackgroundRectangle(velocityrelation)
        velgp = VGroup(velbrec, velocityrelation)
        self.add_foreground_mobject(velgp)
        velgp.move_to(3 * UP)
        self.play(
            AnimationGroup(
                #Animation(Mobject(), run_time = 4.5),
                ApplyMethod(alp.increment_value, 6 * PI, run_time = 9, rate_func = linear),
                Write(velgp),
                lag_ratio = 0.5
            )
        )
        self.wait()
        self.play(
            AnimationGroup(
                ApplyMethod(alp.increment_value, -6 * PI, run_time = 9, rate_func = linear),
                FadeOutAndShiftDown(velgp),
                lag_ratio = 0.5
            )
        )
        self.wait()
        linesgp.clear_updaters()
        path_gp.clear_updaters()
        self.play(FadeIn(spi), FadeOut(path_gp), FadeOut(linesgp), run_time = 5)
        self.wait()
        tgrp = VGroup()
        templine = Line(10 * LEFT, 10 * RIGHT).rotate_about_origin(spiang).move_to(rad_line.get_end()).set_color(PURPLE)
        tangle = Angle(line.get_end(), rad_line.get_end(), templine.get_end(), radius = 0.375)
        tlabel = TexMobject("\\alpha").scale(0.75)
        tlabel.move_to(rad_line.get_end() + 0.6 * (np.cos(spiang / 2) * RIGHT + np.sin(spiang / 2) * UP))
        val = tlabel.get_center() - rad_line.get_end()
        tgrp.add(line, rad_line, templine, tangle, tlabel)
        def tangentupdater(obj):
            trrline, trline, tline, tale, tlbl = obj
            trrline.become(Line(ORIGIN, 10 * RIGHT))
            trrline.rotate_about_origin(omega * alp.get_value())
            trrline.fade(0.75)
            leng = np.exp(omega * alp.get_value() / np.tan(spiang))
            temprli = Line(ORIGIN, leng * RIGHT).set_color(GREEN)
            temprli.rotate_about_origin(omega * alp.get_value())
            trline.become(temprli)
            temprli = Line(10 * LEFT, 10 * RIGHT).rotate_about_origin(spiang).shift(leng * RIGHT).set_color(PURPLE)
            temprli.rotate_about_origin(omega * alp.get_value())
            tline.become(temprli)
            tale.become(Angle(trrline.get_end(), trline.get_end(), tline.get_end(), radius = 0.375))
            testline = Line(ORIGIN, val).shift(leng * RIGHT).rotate_about_origin(alp.get_value())
            tlbl.move_to(testline.get_end())
        self.play(Write(tgrp))
        self.wait()
        tgrp.add_updater(tangentupdater)
        self.add(tgrp)
        self.play(alp.set_value, PI / 4)
        self.wait()
        self.play(alp.increment_value, PI / 3)
        self.wait()
        self.play(alp.increment_value, PI / 2)
        self.play(alp.set_value, 0, run_time = 3)
        self.wait()
        tgrp.clear_updaters()
        self.play(FadeOut(tangle), FadeOut(tlabel))
        self.wait()
        spiang = np.arctan(3.25)
        #self.play(Transform(spi, LogSpiral(start_val = -25, end_val = 25, alpha = spiang)))
        scale_val = 1 / 3
        self.play(
            Transform(spi, LogSpiral(start_val = -25, end_val = 25, alpha = spiang)),
            templine.rotate, -9.602728969052365 * DEGREES, {"about_point": rad_line.get_end()},
            #Rotating(templine, angle = 9.602728969052365 * DEGREES, about_point = rad_line.get_end()),
            self.camera_frame.scale, scale_val,
            self.camera_frame.move_to, rad_line.get_end(),
            run_time = 2
        )
        mag = 1
        linvel = rad_line.get_length() * omega / np.tan(spiang)
        linvec = Arrow(ORIGIN, mag * linvel * RIGHT, max_tip_length_to_length_ratio = 0.25 * scale_val).set_color(YELLOW)
        totvec = Arrow(ORIGIN, mag * rad_line.get_length() * omega * UP, max_tip_length_to_length_ratio = 0.25 * scale_val).set_color(PINK)
        tangvec = Arrow(ORIGIN, mag * rad_line.get_length() * omega * UP + mag * linvel * RIGHT, max_tip_length_to_length_ratio = 0.25 * scale_val).set_color(PURPLE)
        vecgrp = VGroup(linvec, totvec, tangvec)
        #vecgrp.scale(1 / 2, about_point = ORIGIN)
        vecgrp.shift(rad_line.get_length() * RIGHT)
        vecgrp.rotate_about_origin(alp.get_value())
        self.play(Write(tangvec))
        self.wait()
        veclabels = VGroup(
            TexMobject("v_l=k\\cdot r").scale(1 / 3).next_to(linvec.get_end(), direction = UR, buff = 0.05),
            TexMobject("v_n=r\\cdot \\omega").scale(1 / 3).next_to(totvec, direction = LEFT, buff = 0.05),
        )
        rlabel = TexMobject("r").scale(1 / 3).next_to(rad_line, direction = UP, buff = 0.05)
        self.play(Write(rlabel), ApplyWave(rad_line))
        self.wait()
        self.play(Write(linvec), Write(veclabels[0]))
        self.wait()
        self.play(Write(totvec), Write(veclabels[1]))
        self.wait()
        totvec_cpy = totvec.copy().fade(0.5)
        self.play(
            totvec_cpy.shift, linvec.get_length() * RIGHT,
            veclabels[1].next_to, totvec_cpy, {"direction": RIGHT, "buff": linvec.get_length() + 0.05}
        )
        self.wait()
        tang = Angle(linvec.get_end(), rad_line.get_end(), tangvec.get_end())
        tang.scale(1 / 4, about_point = rad_line.get_end())
        ttlabel = tlabel.copy().move_to(rad_line.get_end())
        ttlabel.scale(1 / 2.5, about_point = rad_line.get_end())
        ttlabel.shift(0.225 * (np.cos(spiang / 2) * RIGHT + np.sin(spiang / 2) * UP))
        self.play(Write(tang), Write(ttlabel))
        self.wait()
        angleeqn = TexMobject("\\displaystyle \\tan \\alpha = \\frac {r \\omega}{k r} = \\frac{\\omega}{k} = \\text{const.}").scale(1 / 3)
        angleeqnrect = BackgroundRectangle(angleeqn)
        angleeqngrp = VGroup(angleeqnrect, angleeqn)
        self.add_foreground_mobject(angleeqngrp)
        angleeqngrp.move_to(2 * RIGHT + 0.5 * DOWN)
        self.play(Write(angleeqn))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [angleeqngrp, tang, ttlabel, rlabel, vecgrp, totvec_cpy, veclabels]],
            self.camera_frame.move_to, ORIGIN,
            self.camera_frame.scale, 3,
        )
        self.wait()
        self.play(Rotating(spi, angle = PI, about_point = ORIGIN, rate_func = smooth, run_time = 5))
        self.play(Rotating(spi, angle = -2 * PI, about_point = ORIGIN, rate_func = smooth, run_time = 5))
        self.play(Rotating(spi, angle = PI, about_point = ORIGIN, rate_func = smooth, run_time = 5))
        self.wait()
        spig = LogSpiral(start_val = -50, end_val = 0, alpha = spiang)
        self.add(spig)
        self.play(Indicate(spig, scale_factor = 2), run_time = 2)
        self.play(
            FadeOut(spi),
            self.camera_frame.scale, 1 / 3,
            self.camera_frame.move_to, 0.5 * DOWN,
        )
        rotgrp = VGroup(rad_line, line, templine)
        def rotupdater(obj):
            rl, ll, tl = obj
            ll.become(Line(ORIGIN, 10 * RIGHT))
            ll.rotate_about_origin(omega * alp.get_value())
            ll.fade(0.75)
            leng = np.exp(omega * alp.get_value() / np.tan(spiang))
            temprli = Line(ORIGIN, leng * RIGHT).set_color(GREEN)
            temprli.rotate_about_origin(omega * alp.get_value())
            rl.become(temprli)
            temprli = Line(10 * LEFT, 10 * RIGHT).rotate_about_origin(spiang).shift(leng * RIGHT).set_color(PURPLE)
            temprli.rotate_about_origin(omega * alp.get_value())
            tl.become(temprli)
        rotgrp.add_updater(rotupdater)
        self.add(rotgrp)
        self.play(alp.set_value, -2 * PI, rate_func = there_and_back, run_time = 10)
        self.wait()
        rotgrp.clear_updaters()
        planrotgrp = VGroup(rad_line, line, templine, spig, nplane, polaraxes)
        tspig = spig.copy()
        d = rad_line.get_length()
        def rotoupdater(obj):
            rl, ll, tl, sp, xy, plr = obj
            ll.become(Line(ORIGIN, 10 * RIGHT))
            ll.rotate_about_origin(omega * alp.get_value())
            ll.fade(0.75)
            leng = np.exp(omega * alp.get_value() / np.tan(spiang))
            temprli = Line(ORIGIN, leng * RIGHT).set_color(GREEN)
            temprli.rotate_about_origin(omega * alp.get_value())
            rl.become(temprli)
            temprli = Line(10 * LEFT, 10 * RIGHT).rotate_about_origin(spiang).shift(leng * RIGHT).set_color(PURPLE)
            temprli.rotate_about_origin(omega * alp.get_value())
            tl.become(temprli)
            temprli, tempxy, tempplr = tspig.copy(), tnplane.copy(), tpolaraxes.copy()
            VGroup(rl, ll, tl, temprli, tempxy, tempplr).rotate_about_origin(-omega * alp.get_value())
            VGroup(rl, ll, tl, temprli, tempxy, tempplr).shift((d - leng) * RIGHT)
            for ob, bj in zip([sp, xy, plr], [temprli, tempxy, tempplr]):
                ob.become(bj)
        planrotgrp.add_updater(rotoupdater)
        self.add(planrotgrp)
        self.play(alp.set_value, -PI, rate_func = there_and_back, run_time = 15)
        self.wait()
        planrotgrp.clear_updaters()
        self.add(planrotgrp)
        arcline = Line(rad_line.get_end(), rad_line.get_end() + 0.001 * DOWN).set_color(BLUE)
        theta = -0.5 * PI
        oline = Line(ORIGIN, np.exp(omega * theta / np.tan(spiang)) * RIGHT).set_color(ORANGE)
        oline.rotate_about_origin(theta)
        tangentdot = Dot(arcline.get_end()).scale(1 / 4).set_color(BLUE)
        arcgp = VGroup(spig, arcline, rad_line.copy(), oline, tangentdot)
        self.play(
            FadeOut(nplane), FadeOut(rad_line), FadeIn(tangentdot),
            templine.rotate, -spiang, {"about_point": rad_line.get_end()},
            arcgp.rotate, -spiang, {"about_point": rad_line.get_end()},
            self.camera_frame.move_to, ORIGIN, run_time = 3
        )
        initial_radline = arcgp[2].copy()
        self.add(initial_radline)
        def arcupd(obj):
            sp, al, rl, ol, td = obj
            tempsp = LogSpiral(start_val = -50, end_val = alp.get_value(), alpha = spiang)
            temprl = Line(ORIGIN, tempsp.points[-1]).set_color(GREEN)
            tempal = (d - temprl.get_length()) / np.cos(spiang)
            tempal *= (np.cos(spiang) * LEFT + np.sin(spiang) * DOWN)
            tempal = Line(rad_line.get_end(), rad_line.get_end() + tempal).set_color(BLUE)
            tempol = Line(ORIGIN, np.exp(omega * theta / np.tan(spiang)) * RIGHT).set_color(ORANGE)
            tempol.rotate_about_origin(theta)
            VGroup(tempsp, temprl, tempol).rotate_about_origin(-alp.get_value())
            VGroup(tempsp, temprl, tempol).shift(tempal.get_end() - temprl.get_end())
            VGroup(tempsp, temprl, tempal, tempol).rotate(-spiang, about_point = rad_line.get_end())
            sp.become(tempsp)
            rl.become(temprl)
            al.become(tempal)
            ol.become(tempol)
            td.move_to(al.get_end())
        arcgp.add_updater(arcupd)
        self.add(arcgp)
        self.play(alp.set_value, theta / 2, rate_func = smooth, run_time = 2.5)
        self.wait()
        samp_lines = VGroup()
        perps = VGroup()
        for i in [3700, 3800, 3900]:
            temp = Arrow(arcline.get_end(), spig.points[i], max_tip_length_to_length_ratio = 0.25 * scale_val).set_color(PURPLE)
            samp_lines.add(temp)
        self.play(Write(samp_lines))
        self.wait()
        for i in samp_lines:
            temp = i.copy()
            temp.rotate(PI / 2, about_point = i.get_start())
            #temp.scale(1 / i.get_length(), about_point = i.get_start())
            temp.move_to(i.get_end() + temp.get_center() - temp.get_start())
            temp.become(Arrow(temp.get_start(), temp.get_start() + scale_val * (temp.get_end() - temp.get_start()) / temp.get_length(), max_tip_length_to_length_ratio = 0.25 * scale_val).set_color(RED))
            perps.add(temp)
            self.play(Write(temp))
        self.wait()
        polemov = Arrow(arcgp[2].get_end(), arcgp[2].get_start())
        polemov.rotate(PI / 2, about_point = arcline.get_end())
        polemov.become(Arrow(polemov.get_start(), polemov.get_start() + scale_val * (polemov.get_end() - polemov.get_start()) / polemov.get_length(), max_tip_length_to_length_ratio = 0.25 * scale_val))
        polemov.move_to(arcgp[2].get_start() + polemov.get_center() - polemov.get_start())
        rangle = Angle(arcline.get_end(), polemov.get_start(), polemov.get_end(), radius = 1 / 12).set_color(YELLOW)
        self.play(FadeOut(samp_lines), FadeOut(perps), Write(polemov), Write(rangle),
            #FadeOut(oline),
        )
        self.wait()
        impline = Line(polemov.get_start() + 10 * (polemov.get_end() - polemov.get_start()), polemov.get_start() - 10 * (polemov.get_end() - polemov.get_start())).fade(0.5)
        self.play(Write(impline))
        self.wait()
        self.play(
            FadeOut(rangle),
            FadeOut(polemov),
            #FadeIn(oline),
        )
        self.wait()
        #self.add(arcgp)
        self.play(alp.set_value, theta, rate_func = smooth, run_time = 2.5)
        #self.play(alp.set_value, 3 * theta / 4, rate_func = smooth, run_time = 2.5)
        self.wait()
        horline = Line(arcgp[2].get_start(), arcgp[2].get_start() + 10 * RIGHT)
        #intpt = line_intersection(initial_radline, horline)
        intpt = Intersection(horline, initial_radline)[0].get_center()
        horline.become(Line(arcgp[2].get_start(), intpt).set_color(PURPLE))
        rangle = Angle(initial_radline.get_end(), initial_radline.get_start(), oline.get_start(), radius = 1 / 12)
        tangle = Angle(initial_radline.get_start(), intpt, oline.get_start(), radius = 1 / 12)
        self.play(Write(horline), Write(rangle), Write(tangle))
        self.wait()
        textlabels = VGroup(
            TexMobject("\\alpha"),
            TexMobject("l"),
            TexMobject("r_1"),
            TexMobject("r_0")
        )
        for obj in textlabels:
            obj.scale(1 / 4)
        textlabels[0].next_to(tangle, direction = UL, buff = 0.05 / 4)
        textlabels[1].next_to(horline, direction = DOWN, buff = 0.05 / 2)
        textlabels[2].next_to(oline.get_center(), direction = LEFT, buff = 0.05)
        textlabels[3].next_to(initial_radline.get_center(), direction = RIGHT, buff = 0.05)
        self.play(Write(textlabels))
        arcgp.clear_updaters()
        #scaled_radline = arcgp[2].copy()
        self.remove(arcgp[2])
        arcleneqn = TexMobject("l \\cos(\\alpha)=r_0-r_1").scale(0.5).shift(0.5 * DOWN)
        self.play(Write(arcleneqn))
        self.wait()
        self.play(Transform(arcleneqn, TexMobject("\\displaystyle l=\\frac{r_0-r_1}{\\cos(\\alpha)}").scale(0.5).shift(0.5 * DOWN)))
        self.wait()
        arcgp = VGroup(*[arcgp[i] for i in [0, 1]])
        fadedtspig = tspig.copy().fade(0.5)
        self.play(
            *[FadeOut(obj) for obj in [arcleneqn, textlabels, oline, tangle, rangle, horline, impline, initial_radline]],
            FadeIn(nplane), FadeIn(rad_line), FadeIn(fadedtspig), FadeOut(tangentdot),
            templine.rotate, spiang, {"about_point": rad_line.get_end()},
            arcgp.rotate, spiang, {"about_point": rad_line.get_end()},
            self.camera_frame.scale, 3,
            #self.camera_frame.move_to, 0.5 * DOWN,
            run_time = 2
        )
        self.wait()
        def arcupd1(obj):
            #sp, al, rl, ol, td = obj
            sp, al = obj
            tempsp = LogSpiral(start_val = -50, end_val = alp.get_value(), alpha = spiang)
            temprl = Line(ORIGIN, tempsp.points[-1]).set_color(GREEN)
            tempal = (d - temprl.get_length()) / np.cos(spiang)
            tempal *= (np.cos(spiang) * LEFT + np.sin(spiang) * DOWN)
            tempal = Line(rad_line.get_end(), rad_line.get_end() + tempal).set_color(BLUE)
            #tempol = Line(ORIGIN, np.exp(omega * theta / np.tan(spiang)) * RIGHT).set_color(ORANGE)
            #tempol.rotate_about_origin(theta)
            #VGroup(tempsp, temprl, tempol).rotate_about_origin(-alp.get_value())
            VGroup(tempsp, temprl).rotate_about_origin(-alp.get_value())
            #VGroup(tempsp, temprl, tempol).shift(tempal.get_end() - temprl.get_end())
            tempsp.shift(tempal.get_end() - temprl.get_end())
            sp.become(tempsp)
            #rl.become(temprl)
            al.become(tempal)
            #ol.become(tempol)
            #td.move_to(al.get_end())
        arcgp.add_updater(arcupd1)
        self.add(arcgp)
        self.play(alp.set_value, -12.5, run_time = 5, rate_func = smooth)
        self.wait()
        '''anglegp = VGroup(
            Angle(nplane[0].get_end(), rad_line.get_end(), templine.get_end(), radius = 0.375),
            Angle(ORIGIN, rad_line.get_end(), templine.get_start(), radius = 0.375),
        )'''
        anglegp = Angle(nplane[0].get_end(), rad_line.get_end(), templine.get_end(), radius = 0.375)
        arcleneqn = VGroup(
            TexMobject("\\displaystyle \\text{Arc Length} = \\frac{\\text{Radial length}}{\\cos \\alpha}").move_to(2.5 * UP),
            TexMobject("\\alpha").scale(0.75),
            TexMobject("r").scale(0.75),
            TexMobject("l").scale(0.75),
        )
        arcleneqn[1].move_to(rad_line.get_end() + 0.6 * (np.cos(spiang / 2) * RIGHT + np.sin(spiang / 2) * UP))
        arcleneqn[2].next_to(rad_line, direction = UP, buff = 0.1)
        arcleneqn0 = BackgroundRectangle(arcleneqn[0])
        argrp = VGroup(arcleneqn0, arcleneqn[0])
        self.play(Write(anglegp), Write(arcleneqn[1]), Write(arcleneqn[2]))
        self.wait()
        yintercept = Intersection(Line(5 * UP, 5 * DOWN), templine)
        tanseg = Line(rad_line.get_end(), yintercept.get_center()).set_color(BLUE)
        self.play(
            AnimationGroup(
                #Animation(Mobject(), run_time = 4.5),
                ApplyMethod(alp.set_value, 0, run_time = 10, rate_func = there_and_back),
                Write(argrp),
                lag_ratio = 0.5
            )
        )
        self.wait()
        restofspiral = LogSpiral(start_val = 0, end_val = 25, alpha = spiang)
        self.play(WiggleOutThenIn(tanseg))
        self.wait()
        self.play(FadeIn(restofspiral))
        arcgp.clear_updaters()
        self.wait(5)

class ProblemSnapshot(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        self.wait()
        rlines = lines_grp(a, b, c)
        self.play(
            self.camera_frame.move_to, 2 * UR + RIGHT + 0.5 * UP,
        )
        descr = TextMobject("How do we find the average length of these lines?")
        descr.move_to(3 * RIGHT + 5 * UP)
        self.add(triangle, labels[:3], rlines, descr)
        self.wait()
        #self.play(*[FadeOut(obj) for obj in [tgroup, mlabel, am, altitudes[0], labels[3]]])
        self.wait(5)

class RoseCurve(ParametricFunction):
    CONFIG = {
        "radius"        : 1,
        "start_theta"   : 0,
        "kval"          : 3 / 2,
        "end_theta"     : 4 * PI,
        #"step_size"     : 0.001,
        "density"       : 50 * DEFAULT_POINT_DENSITY_1D,
        "color"         : BLUE,
    }
    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        ParametricFunction.__init__(self, self.pos_func, **kwargs)

    def pos_func(self, t):
        T = self.start_theta + t * (self.end_theta - self.start_theta)
        return self.radius * np.array([
            np.cos(self.kval * T) * np.cos(T),
            np.cos(self.kval * T) * np.sin(T),
            0
        ])

class ProblemIntroOutro(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        c = Circle()
        self.play(Write(c))
        theta = ValueTracker(0)
        sect = Sector(angle = theta.get_value(), color = ORANGE, fill_opacity = 0.125)
        def secupdater(obj):
            obj.become(Sector(angle = theta.get_value(), color = ORANGE, fill_opacity = 0.125))
        sect.add_updater(secupdater)
        self.add(sect)
        self.play(theta.set_value, 2 * PI, run_time = 3, rate_func = smooth)
        sect.clear_updaters()
        self.wait()
        self.play(FadeOut(sect))
        lowerx, higherx = -10, 10
        curveslist = VGroup(
            self.get_graph(lambda t: t, x_min = lowerx, x_max = higherx),
            self.get_graph(lambda t: np.sin(t), x_min = lowerx, x_max = higherx),
            self.get_graph(lambda t: np.cos(t), x_min = lowerx, x_max = higherx),
            RoseCurve(kval = 5 / 2),
            self.get_graph(lambda t: np.exp(t), x_min = lowerx, x_max = higherx),
            Ellipse(),
            self.get_graph(lambda t: np.exp(-t * t / 2), x_min = lowerx, x_max = higherx),
            self.get_graph(lambda t: t ** 2, x_min = lowerx, x_max = higherx),
            Cycloid(start_theta = -4 * PI, end_theta = 4 * PI),
            self.get_graph(lambda t: t ** 3, x_min = lowerx, x_max = higherx),
            self.get_graph(lambda t: t ** (3 / 2), x_min = 0.001, x_max = higherx),
            LogSpiral(start_val = -3 * PI, end_val = 3 * PI),
            self.get_graph(lambda t: np.sinh(t), x_min = lowerx, x_max = higherx),
            self.get_graph(lambda t: 1 / t, x_min = 0.1, x_max = higherx),
            RoseCurve(kval = 7 / 2),
            self.get_graph(lambda t: np.cosh(t), x_min = lowerx, x_max = higherx),
            Cycloid(start_theta = -4 * PI, end_theta = 4 * PI),
        )
        for obj in curveslist:
            self.play(Write(obj))
            self.play(obj.fade, 0.25)
        self.wait()
        curveslist.add(c)
        desctext = TextMobject("Geometry should not include lines (curves) that are like strings, in that they are sometimes straight and sometimes curved, since the ratios between straight and curved lines are not known, and I believe cannot be discovered by human minds, and therefore no conclusion based upon such ratios can be accepted as rigorous and exact.").scale(0.75)
        self.play(
            curveslist.fade, 0.875,
            Write(desctext)
        )
        self.wait()
        self.play(FadeOut(curveslist), FadeOut(desctext))
        self.wait()
        spi = LogSpiral(start_val = -6 * PI, end_val = 6 * PI, alpha = 82.5 * DEGREES).fade(0.875)
        galtext = TextMobject("Who is so blind as not to see that, if there are two equal straight lines, one of which is then bent into a curve, that curve will be equal to the straight line?").scale(0.75)
        self.play(ShowCreation(spi), Write(galtext))
        self.wait(5)


if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "ProblemIntroduction" + " -pl"
    command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 16,17"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    #command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)
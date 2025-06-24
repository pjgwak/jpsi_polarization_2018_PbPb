from JpsiFitter import JpsiFitter

analysis = JpsiFitter(config_module_name='mc_mass') # fitting step
analysis.run()


# = builder = #
# 내부적으로 sigModel, bkgModel 이름을 클래스가 가지고 있는 게 좋을 듯?
# 코드 짜가려면 규칙성 있게 만들어서 나중에 이용자가 코드 짤 때 이해하기 쉽도록 할 것.

# = fitter = #
# # 2. cut 적용하기
# 매개변수 범위 바꾸기 -> 절반만
# Fit option -> 절반만

print('===== Finish runMe.py =====')